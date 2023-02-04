"""
    random_on_sphere()

Draws a random orientation on the unit sphere as a length-3 vector

The z-orientation is drawn a random number between -1 and 1.
The angle in the x-y plane is drawn as a random number between 0 and 2π.
This results in an unbiased random distribution across the sphere.
"""
function random_on_sphere()
    z = rand(Float) * 2 - 1
    r = sqrt(1 - z*z)
    theta = rand(Float) * Float(2 * π)
    (s, c) = sincos(theta)
    return SA[
        r * s,
        r * c,
        z
    ]
end

correct_for_timestep(probability, timestep) = 1 - (1 - probability)^sqrt(timestep)

"""
    random_gauss(diffusivity, timestep)

Draws a random step from a Gaussian distribution with a variance of 2 * diffusivity * timestep.
"""
random_gauss(diffusivity :: Float, timestep :: Float) = sqrt(2 * timestep * diffusivity) * randn(PosVector)

"""
    draw_step!(spin::Spin, diffusivity, timestept, default_properties, [, geometry])

Updates the spin based on a random movement through the given geometry for a given `timestep`:
- draws the next location of the particle after `timestep` with given `diffusivity`.  
  If provided, this displacement will take into account the obstructions in `geometry`.
- the spin orientation might be affected by collisions with the obstructions in `geometry`.
"""
function draw_step!(spin :: Spin, diffusivity :: Float, timestep :: Float, default_properties::GlobalProperties)
    @spin_rng spin begin
        spin.position += random_gauss(diffusivity, timestep)
    end
end
draw_step!(spin :: Spin, diffusivity :: Float, timestep :: Float, default_properties::GlobalProperties, geometry :: Geometry{0}) = draw_step!(spin, diffusivity, timestep, default_properties)
function draw_step!(spin :: Spin{N}, diffusivity :: Float, timestep :: Float, default_properties::GlobalProperties, geometry::Geometry) where {N}
    current_pos = spin.position
    new_pos = random_gauss(diffusivity, timestep) + current_pos
    collision = empty_collision

    @spin_rng spin begin
        for _ in 1:10000
            collision = detect_collision(
                Movement(current_pos, new_pos, one(Float)),
                geometry,
                collision
            )
            if collision === empty_collision
                break
            end

            transfer!.(spin.orientations, correct_for_timestep(MT_fraction(collision, default_properties), timestep))

            permeability_prob = correct_for_timestep(permeability(collision, default_properties), timestep)
            if iszero(permeability_prob) || rand() > permeability_prob
                direction = new_pos .- current_pos
                reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction

                current_pos = collision.distance .* new_pos .+ (1 - collision.distance) .* current_pos
                new_pos = current_pos .+ reflection / norm(reflection) * norm(direction) * (1 - collision.distance)
            else
                current_pos = collision.distance .* new_pos .+ (1 - collision.distance) .* current_pos
                collision = Collision(collision.distance, -collision.normal, collision.properties, collision.index, !collision.inside)
            end
        end
    end
    if collision !== empty_collision
        error("Bounced single particle for 10000 times in single step; terminating!")
    end
    spin.position = new_pos
end

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



