

correct_for_timestep(probability, timestep) = 1 - (1 - probability)^sqrt(timestep)

"""
    random_gauss(diffusivity, timestep)

Draws a random step from a Gaussian distribution with a variance of 2 * diffusivity * timestep.
"""
random_gauss(diffusivity :: Float, timestep :: Float) = sqrt(2 * timestep * diffusivity) * randn(PosVector)

"""
    correct_collision(movement, geometry)

Splits the given movement from point A to point B into multiple steps that bounce off the given obstructions.
This function assumes perfect reflection rather than the diffuse reflection used in [`draw_step`](@ref).
It is used to test the collision detection and resolution, but not actually used in the simulations.
"""
correct_collisions(to_try :: Movement, geometry) = correct_collisions(to_try, Geometry(geometry))
correct_collisions(to_try :: Movement, geometry::Geometry{0}) = [to_try]
function correct_collisions(to_try :: Movement, geometry::Geometry)
    steps = Movement[]
    collision = empty_collision
    while true
        collision = detect_collision(to_try, geometry, collision)
        if collision === empty_collision
            push!(steps, to_try)
            break
        end
        new_pos = collision.distance * to_try.destination + (1 - collision.distance) * to_try.origin
        push!(steps, Movement(
            to_try.origin,
            new_pos,
            collision.distance * to_try.timestep
        ))
        if length(steps) >= 100
            error()
        end
        direction = to_try.destination .- to_try.origin
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction
        new_dest = new_pos .+ reflection / norm(reflection) * norm(direction) * (1 - collision.distance)
        to_try = Movement(
            new_pos,
            new_dest,
            (1 - collision.distance) * to_try.timestep
        )
    end
    steps
end



