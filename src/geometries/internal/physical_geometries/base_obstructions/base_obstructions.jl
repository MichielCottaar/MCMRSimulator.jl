module BaseObstructions

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, get_intersection_params, has_inside, has_single_inside, inside_indices_eltype, isinside_single
import ..PhysicalGeometries: to_property_index
import ...InternalBoundingBoxes: InternalBoundingBox
import ...InternalBoundingBoxes
import ...Properties: GeometryLeafProperties
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh, _mesh_result
import Distributions: Poisson
import Random: rand

abstract type BaseObstruction{N} <: PhysicalGeometry{N} end

child_type(::Type{<:BaseObstruction}) =
    throw(ArgumentError("base obstructions do not have child geometries"))

inside_indices_eltype(::Type{<:BaseObstruction}) = Tuple{}

function surface_sampling end

include("infinite_walls.jl")
include("rounds.jl")
include("overlapping_rounds.jl")
include("triangles.jl")


function random_surface_positions(
    geometry::BaseObstruction{N}, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}, scale_density,
) where {N}
    positions = surface_sampling(geometry, density, scale_density)
    inside = rand(Bool, length(positions))
    return positions, [(i, ) for i in inside]
end

to_property_index(geometry::BaseObstruction, index) = ()

end
