module BaseObstructions

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, child_type, find_intersection, get_intersection_params, has_inside, has_single_inside, inside_indices_eltype, isinside_single, flip
import ..PhysicalGeometries: intersection_index_length, inside_index_length, all_equal_inside_depth
import ...InternalBoundingBoxes: InternalBoundingBox
import ...Indices: ObstructionIndex
import ...InternalBoundingBoxes
import ...Properties: GeometryLeafProperties
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh, _mesh_result
import Distributions: Poisson
import Random: rand

abstract type BaseObstruction{N} <: PhysicalGeometry{N} end

child_type(::Type{<:BaseObstruction}) =
    throw(ArgumentError("base obstructions do not have child geometries"))

inside_indices_eltype(::Type{<:BaseObstruction}) = Tuple{}

intersection_index_length(::Type{<:BaseObstruction}) = 0
inside_index_length(::Type{<:BaseObstruction}) = 0
all_equal_inside_depth(::Type{<:BaseObstruction}) = true

function surface_sampling end

include("infinite_walls.jl")
include("rounds.jl")
include("overlapping_rounds.jl")
include("triangles.jl")


function random_surface_positions(
    geometry::BaseObstruction{N}, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}, scale_density,
) where {N}
    positions, normals = surface_sampling(geometry, density, bounding_box, scale_density)
    intersections = Intersection{N}[]
    for normal in normals
        intersection = Intersection(0.0, normal, true, ObstructionIndex(), false)
        push!(intersections, rand(Bool) ? intersection : flip(intersection))
    end
    positions, intersections
end

end
