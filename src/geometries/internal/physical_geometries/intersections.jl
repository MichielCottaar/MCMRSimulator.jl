module Intersections
import StaticArrays: SVector
import ...Indices: ObstructionIndex, add_index, remove_index

"""
    Intersection{N}(distance, normal, inside, obstruction_index, hit_gap)

Intersection between a path and a physical geometry. `distance` is normalized
to the path from start to destination and `obstruction_index` identifies the
obstruction that was hit.
"""
struct Intersection{N, M}
    distance::Float64
    normal::SVector{N, Float64}
    inside::Bool
    obstruction_index::ObstructionIndex{M}
    hit_gap::Bool
end

Intersection{N}() where {N} = Intersection(Inf, zero(SVector{N, Float64}), false, ObstructionIndex(), false)
function Base.isempty(intersection::Intersection)
    intersection.distance < 0 || intersection.distance > 1
end

add_index(intersection::Intersection, new_index::Int) = Intersection(
    intersection.distance,
    intersection.normal,
    intersection.inside,
    add_index(intersection.obstruction_index, new_index),
    intersection.hit_gap
)
function remove_index(intersection::Intersection)
    index, obstruction_index = remove_index(intersection.obstruction_index)
    return (index, Intersection(
    intersection.distance,
    intersection.normal,
    intersection.inside,
    obstruction_index,
    intersection.hit_gap
))
end

end
