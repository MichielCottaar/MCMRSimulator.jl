module BaseObstructions

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices
import ...InternalBoundingBoxes: InternalBoundingBox
import ...Indices: ObstructionIndex
import ...InternalBoundingBoxes
import ...Properties: GeometryLeafProperties
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh, _mesh_result
import Distributions: Poisson
import Random: rand

abstract type BaseObstruction{N} <: PhysicalGeometry{N} end

function surface_sampling end

include("infinite_walls.jl")
include("rounds.jl")
include("overlapping_rounds.jl")
include("triangles.jl")

function surface_samples_to_intersection(
    positions::AbstractVector{SVector{N, Float64}},
    normals::AbstractVector{SVector{N, Float64}},
) where {N}
    intersections = Intersection{N}[]
    for normal in normals
        inside = rand(Bool)
        push!(intersections, Intersection(0.0, inside ? normal : -normal, inside,
            ObstructionIndex(), false))
    end
    positions, intersections
end

function random_surface_positions(
    geometry::BaseObstruction{N}, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}, scale_density,
) where {N}
    surface_samples_to_intersection(surface_sampling(geometry, density, bounding_box, scale_density)...)
end

size_scale(::InfiniteWall) = Inf
size_scale(round::Round) = round.radius
size_scale(round::OverlappingRound) = round.radius

end
