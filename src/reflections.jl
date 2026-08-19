"""Collision and movement state for diffusing spins."""
module Reflections

import StaticArrays: SVector
import LinearAlgebra: ⋅, norm
import ..Geometries.Internal: Intersection, flip

"""State carried by a spin while it is reflecting from or bound to a surface."""
struct Reflection
    intersection::Union{Nothing, Intersection}
    inside::Bool
    direction::SVector{3, Float64}
    ratio_displaced::Float64
    time_moved::Float64
    distance_moved::Float64
end

function Reflection(
    collision::Intersection,
    direction::SVector{3},
    ratio_displaced,
    time_moved,
    distance_moved,
    permeable=false,
)
    updated_collision = permeable ? flip(collision) : collision
    if permeable
        reflection_direction = direction
        inside = updated_collision.inside
    else
        reflection_direction =
            -2 * (collision.normal ⋅ direction) * collision.normal /
            norm(collision.normal)^2 + direction
        inside = updated_collision.inside
    end

    normalized_direction = reflection_direction / norm(reflection_direction)
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
    nothing,
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
    indices = (geometry_index, obstruction_index)
    collision = Intersection(0., indices, indices, zero(SVector{3, Float64}), inside, false)
    Reflection(collision, inside, direction, ratio_displaced, time_moved, distance_moved)
end

has_intersection(reflection::Reflection) = !isnothing(reflection.intersection)
has_intersection(intersection::Intersection) = true

has_hit(reflection::Reflection) = isnothing(reflection.intersection) ? () : reflection.intersection.indices
previous_hit(reflection::Reflection) = isnothing(reflection.intersection) ? nothing :
    (reflection.intersection.indices..., reflection.intersection.inside, reflection.intersection.distance)

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
