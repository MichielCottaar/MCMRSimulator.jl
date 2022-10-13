"Intermediate object used internally to represent a movement from one position to another"
struct Movement{N, M}
    origin :: SVector{N, Float}
    destination :: SVector{N, Float}
    orientations :: SVector{M, SpinOrientation}
    timestep :: Float
end

function Movement(origin::AbstractArray, destination::AbstractArray, timestep::Real) 
    ndim = length(origin)
    Movement{ndim, 0}(SVector{ndim, Float}(origin), SVector{ndim, Float}(destination), SVector{0, SpinOrientation}(), Float(timestep))
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
    random_gauss(diffusivity, timestep)

Draws a random step from a Gaussian distribution with a variance of 2 * diffusivity * timestep.
"""
random_gauss(diffusivity :: Float, timestep :: Float) = sqrt(2 * timestep * diffusivity) * randn(PosVector)

"""
    draw_step(current::Spin, diffusivity::Float, timestep::Float[, geometry::Obstructions])

Returns a new spin with the state of the particle after `timestep`:
- draws the next location of the particle after `timestep` with given `diffusivity`.  
  If provided, this displacement will take into account the obstructions in `geometry`.
- the spin orientation might be affected by collisions with the obstructions in `geometry`.
"""
function draw_step(current :: Spin{N}, diffusivity :: Float, timestep :: Float) where {N}
    new_rng = @spin_rng current begin
        new_pos = current.position .+ random_gauss(diffusivity, timestep)
    end
    Spin{N}(new_pos, current.orientations, new_rng)
end
draw_step(current :: Spin, diffusivity :: Float, timestep :: Float, geometry :: Tuple{}) = draw_step(current, diffusivity, timestep)
draw_step(current :: Spin, diffusivity :: Float, timestep :: Float, geometry :: AbstractVector{<:Obstruction}) = draw_step(current, diffusivity, timestep, Tuple(geometry))
function draw_step(current :: Spin{N}, diffusivity :: Float, timestep :: Float, geometry::Tuple) where {N}
    proposed = draw_step(current, diffusivity, timestep)
    displacement = norm(current.position .- proposed.position)
    empty = empty_collision(N)
    collision = empty

    current_pos = current.position
    new_pos = proposed.position
    orient = proposed.orientations

    final_rng = @spin_rng proposed begin
        for _ in 1:1000
            collision = detect_collision(
                Movement(current_pos, new_pos, orient, one(Float)),
                geometry,
                collision
            )
            if collision === empty
                break
            end
            orient = collision.new_orient
            current_pos = collision.distance .* new_pos .+ (1 - collision.distance) .* current_pos
            direction = random_on_sphere()
            displacement = (1 - collision.distance) * displacement
            new_pos = current_pos .+ (sign(direction ⋅ collision.normal) * displacement) .* direction
        end
    end
    if collision !== empty
        error("Bounced single particle for 1000 times in single step; terminating!")
    end
    Spin(new_pos, orient, final_rng)
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

function correct_collisions(to_try :: Movement{N, M}, geometry::Tuple) where {N, M}
    steps = Movement[]
    empty = empty_collision(M)
    collision = empty
    while true
        collision = detect_collision(to_try, geometry, collision)
        if collision === empty
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
    Collision(distance, normal, new_orient, id[, index])

A detected collision along the movement.

# Parameters
- `distance`: number between 0 and 1 indicating the distance of the collision from the origin (0) to the destination point (1)
- `normal`: normal of the obstruction at the collision site. To get correct reflection the normal should point in the direction of the incoming particle.
- `new_orient`: new spin orientation after magnetisation transfer with the obstruction.
"""
struct Collision{N}
    distance :: Float
    normal :: PosVector
    new_orient :: SVector{N, SpinOrientation}
    id :: UUID
    index :: Int
    Collision(distance, normal, new_orient, id, index) = new{length(new_orient)}(prevfloat(distance), normal, new_orient, id, index)
end

Collision(distance, normal, new_orient, id; index=0) = Collision(distance, normal, new_orient, id, index)


const empty_collision_id = uuid1()
empty_collision(N::Int) = Collision(Inf, SA[0, 0, 0], SVector{N, SpinOrientation}([SpinOrientation(0., 0., 0.) for _ in 1:N]), empty_collision_id, 0)

"""
    detect_collision(movement, obstructions[, previous])

Returns a [`Collision`](@ref) object if the given `movement` crosses any obstructions.
The first collision is always returned.
If no collision is detected, `empty_collision` will be returned
"""
detect_collision(movement :: Movement{M, N}, obstructions) where {M, N} = detect_collision(movement, obstructions, empty_collision(N))
detect_collision(movement :: Movement{M, N}, obstructions :: Tuple{}, previous::Collision{N}) where {M, N} = empty_collision(N)
detect_collision(movement :: Movement, obstructions :: Tuple{<:Obstruction}, previous::Collision) = detect_collision(movement, obstructions[1], previous)

function detect_collision(movement :: Movement{M, N}, obstructions :: Tuple, previous::Collision{N}) where {M, N}
    collision = empty_collision(N)
    for o in obstructions
        c_new = detect_collision(movement, o, previous)
        if c_new.distance < collision.distance
            collision = c_new
        end
    end
    collision
end

