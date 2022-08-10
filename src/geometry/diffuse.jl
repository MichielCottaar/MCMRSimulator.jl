"Intermediate object used internally to represent a movement from one position to another"
struct Movement
    origin :: PosVector
    destination :: PosVector
    timestep :: Float
end

"""
    random_on_sphere()

Draws a random orientation on the unit sphere as a length-3 vector

The z-orientation is drawn a random number between -1 and 1.
The angle in the x-y plane is drawn as a random number between 0 and 2π.
This results in an unbiased random distribution across the sphere.
"""
function random_on_sphere()
    z = rand() * 2. - 1.
    r = sqrt(1. - z*z)
    theta = rand() * 2 * π
    (s, c) = sincos(theta)
    return SA[
        r * s,
        r * c,
        z
    ]
end

"""
    draw_step([current_pos::PosVector], diffusivity::Float, timestep::Float[, geometry::Obstructions])

Draws the next location of the particle after `timestep` with given `diffusivity`.
If provided, this displacement will take into account the obstructions in `geometry`.
If the `current_pos` is not provided only the displacement is returned (will only work for empty `geometry`).
"""
draw_step(diffusivity :: Float, timestep :: Float) = sqrt(2. * timestep * diffusivity) * @SVector randn(3)
draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float) = current_pos .+ draw_step(diffusivity, timestep)
draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float, geometry :: Obstructions{0}) = draw_step(current_pos, diffusivity, timestep)
function draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float, geometry :: Obstructions)
    original_position = current_pos
    new_pos = draw_step(current_pos, diffusivity, timestep)
    displacement = norm(new_pos .- current_pos)
    for _ in 1:1000
        collision = detect_collision(
            Movement(current_pos, new_pos, 1.),
            geometry,
            original_position
        )
        if isnothing(collision)
            return new_pos
        end
        current_pos = collision.distance .* new_pos .+ (1 - collision.distance) .* current_pos
        direction = random_on_sphere()
        displacement = (1 - collision.distance) * displacement
        new_pos = current_pos .+ (sign(direction ⋅ collision.normal) * displacement) .* direction
    end
    error("Bounced single particle for 1000 times in single step; terminating!")
end

"""
    correct_collision(movement, geometry)

Splits the given movement from point A to point B into multiple steps that bounce off the given obstructions.
This function assumes perfect reflection rather than the diffuse reflection used in [`draw_step`](@ref).
It is used to test the collision detection and resolution, but not actually used in the simulations.
"""
correct_collisions(to_try :: Movement, geometry :: Obstruction) = correct_collisions(to_try, SVector{1}(geometry))
correct_collisions(to_try :: Movement, geometry :: Obstructions{0}) = [to_try]
correct_collisions(to_try :: Movement, geometry :: AbstractVector{<:Obstruction}) = correct_collisions(to_try, SVector{length(geometry)}(geometry))

function correct_collisions(to_try :: Movement, geometry :: Obstructions)
    steps = Movement[]
    start = to_try.origin
    while true
        collision = detect_collision(to_try, geometry, start)
        if isnothing(collision) || collision.distance > 1
            push!(steps, to_try)
            break
        end
        new_pos = collision.distance * to_try.destination + (1 - collision.distance) * to_try.origin
        push!(steps, Movement(
            to_try.origin,
            new_pos,
            collision.distance * to_try.timestep
        ))
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


"""
    Collision(distance::Float, normal::PosVector)

A detected collision along the movement.

# Parameters
- `distance`: number between 0 and 1 indicating the distance of the collision from the origin (0) to the destination point (1)
- `normal`: normal of the obstruction at the collision site. To get correct reflection the normal should point in the direction of the incoming particle.
"""
struct Collision
    distance :: Float
    normal :: PosVector
    Collision(distance, normal) = new(distance * (1. - eps(Float)^0.75), normal)
end


"""
    detect_collision(movement, obstructions)

Returns a [`Collision`](@ref) object if the given `movement` crosses any obstructions.
The first collision is always returned.
If no collision is detected, `nothing` will be returned
"""
detect_collision(movement :: Movement, obstructions :: Obstructions{0}, origin::PosVector) = nothing
detect_collision(movement :: Movement, obstructions :: Obstructions{1}, origin::PosVector) = detect_collision(movement, obstructions[1], origin)

function detect_collision(movement :: Movement, obstructions :: Obstructions, origin::PosVector)
    collision = nothing
    for o in obstructions
        c_new = detect_collision(movement, o, origin)
        if !isnothing(c_new) && (isnothing(collision) || c_new.distance < collision.distance)
            collision = c_new
        end
    end
    collision
end

