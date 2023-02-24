"""
    Reflection(collision, movement, timestep)
    Reflection(collision, timestep, fraction_used, surface_pos, direction)

Represents a reflection of a particle after colliding (see [`Collision`](@ref)).
"""
struct Reflection
    collision::Collision
    timestep::Float
    fraction_used::Float
    surface_pos::PosVector
    direction::PosVector
end

function Reflection(collision::Collision, movement::Movement, timestep, fraction_used, permeable)
    surface_pos = @. (movement.origin * (1 - collision.distance) + movement.destination * collision.distance)

    direction = movement.destination - movement.origin
    if permeable
        flipped = Collision(collision.distance, -collision.normal, collision.properties, collision.index, !collision.inside)
        return Reflection(flipped, timestep, fraction_used, surface_pos, direction * (1 - collision.distance))
    else
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction
        new_direction = reflection * (norm(direction) * (1 - collision.distance) / norm(reflection))
        return Reflection(collision, timestep, fraction_used, surface_pos, new_direction)
    end
end

"""
    final_position(reflection[, remaining_time])

Returns the position that the spin should end up if it has `remaining_time` ms left in this timestep.

If `remaining_time` is not provided, the time remaining is assumed to be the time remaining in this timestep.
"""
final_position(reflection::Reflection) = reflection.direction + reflection.surface_pos
function final_position(reflection::Reflection, new_time)
    additional_movement = sqrt(reflection.fraction_used + new_time / reflection.timestep) - reflection.fraction_used
    return reflection.surface_pos + reflection.direction * (additional_movement / (1 - reflection.fraction_used))
end


const empty_reflection = Reflection(empty_collision, zero(Float), zero(Float), zero(PosVector), zero(PosVector))