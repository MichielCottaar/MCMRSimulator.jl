"""Collision and movement state for diffusing spins."""
module Reflections

import StaticArrays: SVector
import LinearAlgebra: ⋅, norm
import ..Geometries.Internal: Intersection, ObstructionIndex

"""State carried by a spin while it is reflecting from or bound to a surface."""
struct Reflection
    intersection::Intersection{3}
    inside::Bool
    direction::SVector{3, Float64}
    ratio_displaced::Float64
    time_moved::Float64
    distance_moved::Float64
end

function Reflection(
    collision::Intersection{3},
    direction::SVector{3},
    ratio_displaced,
    time_moved,
    distance_moved,
    permeable=false,
)
    if permeable
        reflection_direction = direction
        inside = !collision.inside
    else
        reflection_direction =
            -2 * (collision.normal ⋅ direction) * collision.normal /
            norm(collision.normal)^2 + direction
        inside = collision.inside
    end

    normalized_direction = reflection_direction / norm(reflection_direction)
    updated_collision = Intersection(
        collision.distance,
        collision.normal,
        inside,
        collision.obstruction_index,
        collision.hit_gap,
    )
    Reflection(
        updated_collision,
        inside,
        normalized_direction,
        ratio_displaced,
        time_moved,
        distance_moved,
    )
end

"""Create movement state for a free spin at the start of a timestep."""
Reflection(ratio_displaced) = Reflection(
    Intersection{3}(),
    false,
    zero(SVector{3, Float64}),
    ratio_displaced,
    0.,
    0.,
)

"""Temporary constructor for surface samplers that still return flat indices."""
function Reflection(
    geometry_index::Integer,
    obstruction_index::Integer,
    inside::Bool,
    direction::SVector{3},
    ratio_displaced,
    time_moved,
    distance_moved,
)
    index = ObstructionIndex(SVector{2, Int}(geometry_index, obstruction_index))
    collision = Intersection(0., zero(SVector{3, Float64}), inside, index, false)
    Reflection(collision, inside, direction, ratio_displaced, time_moved, distance_moved)
end

has_intersection(reflection::Reflection) = !isempty(reflection.intersection.obstruction_index.indices)
has_intersection(intersection::Intersection) = !Base.isempty(intersection)

has_hit(reflection::Reflection) = reflection.intersection.obstruction_index
previous_hit(reflection::Reflection) = reflection.intersection

function direction(reflection::Reflection, new_time, diffusivity)
    displacement_size =
        reflection.ratio_displaced *
        sqrt(2 * diffusivity * (reflection.time_moved + new_time)) -
        reflection.distance_moved
    @assert displacement_size > 0
    reflection.direction * displacement_size
end

const empty_reflection = Reflection(0.)

end
