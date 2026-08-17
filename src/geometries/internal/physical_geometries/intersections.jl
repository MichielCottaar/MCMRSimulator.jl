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

Intersection{N, M}() where {N, M} = Intersection(Inf, zero(SVector{N, Float64}), false, ObstructionIndex{M}(), false)
Intersection{N}() where {N} = Intersection{N, 0}()
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

"""
    remove_expected_index(intersection, expected_index)

If the index in `intersection` matches `expected_index` an intersection is returned with that element removed.  Otherwise, an empty intersection is returned.
"""
function remove_expected_index(
    intersection::Intersection{N, M},
    expected::Int,
) where {N, M}
    empty_intersection = Intersection{N, M - 1}()
    Base.isempty(intersection) && return empty_intersection
    isempty(intersection.obstruction_index.indices) &&
        throw(ArgumentError("Intersection indices are unexpectedly empty."))
    intersection.obstruction_index.indices[1] == expected || return empty_intersection
    remove_index(intersection)[2]
end

function remove_expected_index(intersection::Intersection{N, 0}, expected::Int) where {N}
    Base.isempty(intersection) && return intersection
    throw(ArgumentError("Intersection indices are unexpectedly empty."))
end

"""
    flip(intersection)

Flips a collision of the spin with a surface to the other side.

It returns a new collision. This should be called if a spin passes through to the other side of the surface due to permeability.
"""
flip(intersection::Intersection) = Intersection(
    intersection.distance,
    -intersection.normal,
    !intersection.inside,
    intersection.obstruction_index,
    intersection.hit_gap
)

end
