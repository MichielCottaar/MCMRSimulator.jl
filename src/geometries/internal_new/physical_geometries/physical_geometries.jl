module PhysicalGeometries
import StaticArrays: SVector
import ..Indices: ObstructionIndex

"""
    PhysicalGeometry{N}

Represents a physical geometry in N-dimensional space. This is an abstract type that can be extended to define specific geometries and their properties.
"""
abstract type PhysicalGeometry{N} end

"""
    Intersection{N}(distance, normal, inside, obstruction_index, hit_gap)

Intersection between a path and a physical geometry. `distance` is normalized
to the path from start to destination and `obstruction_index` identifies the
obstruction that was hit.
"""
struct Intersection{N}
    distance::Float64
    normal::SVector{N, Float64}
    inside::Bool
    obstruction_index::ObstructionIndex
    hit_gap::Bool
end

Intersection{N}() where {N} = Intersection(Inf, zero(SVector{N, Float64}), false, ObstructionIndex(), false)
function Base.isempty(intersection::Intersection)
    intersection.distance < 0 || intersection.distance > 1
end

"""Find the first intersection of a path with `geometry`."""
function detect_intersection end

"""Return whether a geometry has an inside region."""
function has_inside end

"""Return the obstruction indices containing a position."""
function inside_indices end

include("groups.jl")
include("transformations.jl")
include("base_obstructions/base_obstructions.jl")

end
