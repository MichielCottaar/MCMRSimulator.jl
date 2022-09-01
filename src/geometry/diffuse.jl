"Intermediate object used internally to represent a movement from one position to another"
struct Movement{N}
    origin :: SVector{N, Float}
    destination :: SVector{N, Float}
    timestep :: Float
end

function Movement(origin::AbstractArray, destination::AbstractArray, timestep::Real) 
    ndim = length(origin)
    Movement{ndim}(SVector{ndim, Float}(origin), SVector{ndim, Float}(destination), Float(timestep))
end

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

"""
    draw_step([current_pos::PosVector], diffusivity::Float, timestep::Float[, geometry::Obstructions])

Draws the next location of the particle after `timestep` with given `diffusivity`.
If provided, this displacement will take into account the obstructions in `geometry`.
If the `current_pos` is not provided only the displacement is returned (will only work for empty `geometry`).
"""
draw_step(diffusivity :: Float, timestep :: Float) = sqrt(2 * timestep * diffusivity) * randn(PosVector)
draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float) = current_pos .+ draw_step(diffusivity, timestep)
draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float, geometry :: Tuple{}) = draw_step(current_pos, diffusivity, timestep)
draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float, geometry :: AbstractVector{<:Obstruction}) = draw_step(current_pos, diffusivity, timestep, Tuple(geometry))
function draw_step(current_pos :: PosVector, diffusivity :: Float, timestep :: Float, geometry::Tuple)
    new_pos = draw_step(current_pos, diffusivity, timestep)
    displacement = norm(new_pos .- current_pos)
    collision = empty_collision
    for _ in 1:1000
        collision = detect_collision(
            Movement(current_pos, new_pos, one(Float)),
            geometry,
            collision
        )
        if collision === empty_collision
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
correct_collisions(to_try :: Movement, geometry :: Obstruction) = correct_collisions(to_try, tuple(geometry))
correct_collisions(to_try :: Movement, geometry :: Tuple{}) = [to_try]
correct_collisions(to_try :: Movement, geometry :: AbstractVector{<:Obstruction}) = correct_collisions(to_try, tuple(geometry...))

function correct_collisions(to_try :: Movement, geometry::Tuple)
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
    id :: UUID
    index :: Int
    Collision(distance, normal, id, index) = new(prevfloat(distance), normal, id, index)
end

Collision(distance, normal, id; index=0) = Collision(distance, normal, id, index)


const empty_collision = Collision(Inf, SA[0, 0, 0], uuid1(), 0)

"""
    detect_collision(movement, obstructions)

Returns a [`Collision`](@ref) object if the given `movement` crosses any obstructions.
The first collision is always returned.
If no collision is detected, `empty_collision` will be returned
"""
detect_collision(movement :: Movement, obstructions :: Tuple{}, previous=empty_collision::Collision) = empty_collision
detect_collision(movement :: Movement, obstructions :: Tuple{<:Obstruction}, previous=empty_collision::Collision) = detect_collision(movement, obstructions[1], previous)

function detect_collision(movement :: Movement, obstructions :: Tuple, previous=empty_collision::Collision)
    collision = empty_collision
    for o in obstructions
        c_new = detect_collision(movement, o, previous)
        if c_new.distance < collision.distance
            collision = c_new
        end
    end
    collision
end

