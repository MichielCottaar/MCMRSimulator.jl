"""
    Reflection(collision, movement, timestep)
    Reflection(collision, timestep, fraction_used, direction)

Represents a reflection of a particle after colliding (see [`Collision`](@ref)).
"""
struct Reflection
    collision::Collision
    ratio_displaced::Float
    time_moved::Float
    distance_moved::Float
    direction::PosVector
end

function Reflection(collision::Collision, movement::Movement, ratio_displaced, time_moved, distance_moved, permeable)
    direction = movement.destination - movement.origin
    if permeable
        flipped = Collision(collision.distance, -collision.normal, collision.properties, collision.index, !collision.inside)
        return Reflection(flipped, ratio_displaced, time_moved, distance_moved, direction / norm(direction))
    else
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction
        return Reflection(collision, ratio_displaced, time_moved, distance_moved, reflection / norm(reflection))
    end
end

"""
    direction(reflection, remaining_time)

Returns the displacement that the spin should end up if it has `remaining_time` ms left in this timestep.
"""
function direction(reflection::Reflection, new_time)
    displacement_size = reflection.ratio_displaced * sqrt(reflection.time_moved + new_time) - reflection.distance_moved
    @assert displacement_size > 0
    return reflection.direction * displacement_size
end


"""
    new_reflection(ratio_displaced[, collision=empty_collision])

Produces a new `Reflection` with no previous movement.
If no `collision` is provided it is also assumed to have no previous collisions (i.e., a free particle at the start of any timestep).
If a `collision` is provided the particle is supposed to have been released from this surface (i.e., a stuck particle at the start of the simulation).
"""
new_reflection(ratio_displaced, collision=empty_collision) = Reflection(collision, Float(ratio_displaced), zero(Float), zero(Float), zero(PosVector))