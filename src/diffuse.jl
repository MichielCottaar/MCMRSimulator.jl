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
    return SA[
        r * sin(theta),
        r * cos(theta),
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
    new_pos = draw_step(current_pos, diffusivity, timestep)
    displacement = norm(new_pos .- current_pos)
    while true
        collision = detect_collision(
            Movement(current_pos, new_pos, 1.),
            geometry
        )
        if isnothing(collision)
            return new_pos
        end
        flip_normal = (new_pos .- current_pos) ⋅ collision.normal
        current_pos = collision.distance .* new_pos .+ (1 - collision.distance) .* current_pos
        direction = random_on_sphere()
        displacement = (1 - collision.distance) * displacement
        new_pos = current_pos .+ (((flip_normal * (direction ⋅ collision.normal)) > 0 ? -1 : 1) * displacement) .* direction
    end
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
    while true
        collision = detect_collision(to_try, geometry)
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
A detected collision along the movement.

# Parameters
- `distance`: number between 0 and 1 indicating the distance of the collision from the origin (0) to the destination point (1)
- `normal`: normal of the obstruction at the collision site.
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
detect_collision(movement :: Movement, obstructions :: Obstructions{0}) = nothing

function detect_collision(movement :: Movement, obstructions :: Obstructions)
    collision = nothing
    for o in obstructions
        c_new = detect_collision(movement, o)
        if !isnothing(c_new) && (isnothing(collision) || c_new.distance < collision.distance)
            collision = c_new
        end
    end
    collision
end


"""
    Repeated(obstructions, repeats)

Underlying `obstructions` are repeated ad infinitum as described in `repeats` (a length-3 vector).
Zeros or infinity in `repeats` indicate that no repeating should be applied in this direction.
To repeat in a different direction that the cardinal directions first apply `Repeated` and then [`Transformed`](@ref).

If the obstructions are larger than the repeat size they will be cut off!
"""
struct Repeated{O <: Obstructions} <: Obstruction
    obstructions :: O
    repeats :: PosVector
    function Repeated(obstructions, repeats)
        o = isa(obstructions, Obstruction) ? SVector{1}([obstructions]) : SVector{length(obstructions)}(obstructions)
        rs = SVector{3, Float}([iszero(r) ? Inf : abs(r) for r in repeats])
        any(isfinite.(rs)) ? new{typeof(o)}(o, rs) : obstructions
    end
end


function detect_collision(movement :: Movement, repeat :: Repeated)
    origin = movement.origin ./ repeat.repeats .+ 0.5
    destination = movement.destination ./ repeat.repeats .+ 0.5
    for (_, t1, p1, t2, p2) in ray_grid_intersections(origin, destination)
        f(r, p, p_orig) = isfinite(r) ? r * (p - 0.5) : p_orig
        pos1 = f.(repeat.repeats, p1, (movement.destination .* t1) .+ (movement.origin .* (1 - t1)))
        pos2 = f.(repeat.repeats, p2, (movement.destination .* t2) .+ (movement.origin .* (1 - t2)))
        c = detect_collision(
            Movement(pos1, pos2, 1.),
            repeat.obstructions
        )
        if !isnothing(c)
            return Collision(
                c.distance * (t2 - t1) + t1,
                c.normal
            )
        end
    end
    return nothing
end

struct RayGridIntersections
    origin :: PosVector
    destination :: PosVector
    direction :: PosVector
    all_next_hits :: MVector{3, Float}
    current_voxel :: MVector{3, Int}
end

"""
    ray_grid_intersections(origin, destination)

Computes all voxels crossed by a ray between `origin` and `destination` with a 1x1x1 grid.
Both origin and destination are length-3 vectors.
The returned object is an iterator returning a tuple with:
- 3-length vector with the voxel that we are crossing through
- Float with the time the ray entered voxel (0=`origin`, 1=`destination`)
- 3-length vector with position within voxel that the ray entered (i.e., numbers between 0 and 1)
- Float with the time the ray left the voxel (0=`origin`, 1=`destination`)
- 3-length vector with position within voxel that the ray left (i.e., numbers between 0 and 1)
"""
function ray_grid_intersections(origin :: PosVector, destination :: PosVector)
    direction = destination .- origin
    within_voxel = mod.(origin, 1)
    all_next_hits = MVector{3, Float}(map((d, w) -> (d > 0 ? 1 - w : w) / abs(d), direction, within_voxel))
    current_voxel = MVector{3, Int}(map(o->Int(floor(o)), origin))
    return RayGridIntersections(origin, destination, direction, all_next_hits, current_voxel)
end

Base.iterate(rgi::RayGridIntersections) = Base.iterate(rgi, (rgi.origin, zero(Float)))
function Base.iterate(rgi::RayGridIntersections, state::Tuple{PosVector, Float})
    (prev_pos::PosVector, prev_time::Float) = state
    if prev_time >= 1.
        return nothing
    end
    dimension = argmin(rgi.all_next_hits)
    time_to_hit = rgi.all_next_hits[dimension]
    next_time = prev_time + time_to_hit
    if next_time > 1.
        return (
            (PosVector(rgi.current_voxel), prev_time, prev_pos - rgi.current_voxel, Float(1.), rgi.destination - rgi.current_voxel), 
            (rgi.destination, next_time)
        )
    end
    next_pos = rgi.origin .+ (rgi.direction .* next_time)
    res = (
        (PosVector(rgi.current_voxel), prev_time, prev_pos - rgi.current_voxel, next_time, next_pos - rgi.current_voxel),
        (next_pos, next_time)
    )
    for dim in 1:3
        if dim == dimension
            rgi.all_next_hits[dim] = 1 / abs(rgi.direction[dim])
            rgi.current_voxel[dim] += sign(rgi.direction[dim])
        else
            rgi.all_next_hits[dim] -= time_to_hit
        end
    end
    return res
end

Base.length(rgi::RayGridIntersections) = sum(abs.(Int.(floor.(rgi.destination)) .- Int.(floor.(rgi.origin)))) + 1
Base.eltype(::RayGridIntersections) = Tuple{PosVector, Float, PosVector, Float, PosVector}


"""
    Transformed(obstructions::Obstructions, transform::CoordinateTransformations.Transformation)

Underlying `obstructions` are linearly transformed (e.g., rotated or shifted) using the given `transform`.

"""
struct Transformed{O <: Obstructions, T <: CoordinateTransformations.Transformation} <: Obstruction
    obstrucations :: O
    transform :: T
    inverse :: T
    function Transformed(obstrucations, transform :: CoordinateTransformations.Transformation)
        if isa(obstrucations, Transformed)
            return Transformed(obstrucations.obstrucations, transform ∘ obstrucations.transform )
        else
            o = isa(obstrucations, Obstruction) ? SVector{1}([obstrucations]) : SVector{length(obstrucations)}(obstrucations)
            return new{typeof(o), typeof(transform)}(o, transform, inv(transform))
        end
    end
end

function detect_collision(movement :: Movement, transform :: Transformed) 
    c = detect_collision(
        Movement(
            transform.inverse(movement.origin),
            transform.inverse(movement.destination),
            movement.timestep
        ),
        transform.obstrucations
    )
    if isnothing(c)
        return c
    end
    Collision(
        c.distance,
        transform.transform(c.normal) .- transform.transform(zero(PosVector))
    )
end

"""
    Wall([[normal,] offset])

Infinitely large wall with an orientation given by `normal` (default in the x-direction).
The normal can also be defined using :x, :y, or :z, to point in that cardinal direction.
The offset of the wall from the origin along this `normal` is given by `offset` (so that `offset .* normal` is on the wall).
"""
struct Wall <: Obstruction
end

Wall(offset :: Float) = offset == 0 ? Wall() : Transformed(Wall(), CoordinateTransformations.Translation(offset, 0., 0.))


function Wall(normal :: PosVector, offset :: Float)
    n = normal ./ norm(normal)
    shifted = Wall(offset * norm(normal))
    if isapprox(n, SA[1, 0, 0], atol=1e-10)
        return shifted
    end
    rot_axis = cross(SA[1, 0, 0], n)
    rot_angle = acos(n[1])
    Transformed(shifted, CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Wall(sym :: Symbol, offset :: Float)
    direction = Dict(
        :x => SA[1., 0., 0.],
        :y => SA[0., 1., 0.],
        :z => SA[0., 0., 1.],
    )
    Wall(direction[sym], offset)
end

function detect_collision(movement :: Movement, wall :: Wall)
    origin = movement.origin[1]
    destination = movement.destination[1]
    if origin * destination >= 0
        return nothing
    end
    total_length = abs(origin - destination)
    Collision(
        abs(origin) / total_length,
        SA[1, 0, 0]
    )
end

"""
    Sphere([radius[, location]])

Creates a hollow sphere with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
"""
struct Sphere <: Obstruction
    radius :: Float
end

function Sphere(radius :: Real, location :: AbstractVector)
    radius = convert(Float, radius)
    location = convert(SVector{3, Float}, location)
    if all(iszero.(location))
        return Sphere(radius)
    else
        return Transformed(Sphere(radius), CoordinateTransformations.Translation(location))
    end
end

detect_collision(movement :: Movement, sphere :: Sphere) = sphere_collision(movement.origin, movement.destination, sphere.radius)

function sphere_collision(origin :: PosVector, destination :: PosVector, radius :: Float)
    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    if (
        (origin[1] > radius && destination[1] > radius) ||
        (origin[1] < -radius && destination[1] < -radius) ||
        (origin[2] > radius && destination[2] > radius) ||
        (origin[2] < -radius && destination[2] < -radius) ||
        (origin[3] > radius && destination[3] > radius) ||
        (origin[3] < -radius && destination[3] < -radius)
    )
        return nothing
    end
    diff = destination - origin
    a = sum(diff .* diff)
    b = sum(2 .* origin .* diff)
    c = sum(origin .* origin)
    determinant = b ^ 2 - 4 * a * (c - radius ^ 2)
    if determinant < 0
        return nothing
    end
    sd = sqrt(determinant)
    ai = inv(a)

    # first try solution with lower s
    solution = -(b + sd) * 0.5 * ai
    if solution > 1
        return nothing
    elseif solution <= 0
        solution = (-b + sd) * 0.5 * ai
        if solution > 1 || solution < 0
            return nothing
        end
    end
    point_hit = solution * destination + (1 - solution) * origin
    Collision(
        solution,
        point_hit
    )
end

"""
    Cylinder([radius[, orientation[, location]]])

Creates a hollow cylinder with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
The orientation of the cylinder (default: z-direction) can be given as a symbol of the cardinal orientation (:x, :y, or :z) or as a length-3 vector.
"""
struct Cylinder <: Obstruction
    radius :: Float
end

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real})
    radius = Float(radius)
    o = SVector{3, Float}(orientation / norm(orientation))
    if isapprox(o, SA[0, 0, 1], atol=1e-10)
        return Cylinder(radius)
    end
    rot_axis = cross(SA[0, 0, 1], o)
    rot_angle = acos(o[end])
    Transformed(Cylinder(radius), CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real}, location :: AbstractVector{<:Real})
    c = Cylinder(radius, orientation)
    location = SVector{3, Float}(location)
    return all(iszero.(location)) ? c : Transformed(c, CoordinateTransformations.Translation(location...))
end

function Cylinder(radius :: Real, sym :: Symbol, offset :: AbstractVector{<:Real})
    orientation = Dict(
        :x => SA[1., 0., 0.],
        :y => SA[0., 1., 0.],
        :z => SA[0., 0., 1.],
    )
    Cylinder(radius, orientation[sym], offset)
end

function detect_collision(movement :: Movement, cylinder :: Cylinder)
    o = SA[movement.origin[1], movement.origin[2], 0.]
    d = SA[movement.destination[1], movement.destination[2], 0.]
    sphere_collision(o, d, cylinder.radius)
end


"""
    cylinder_plane(radius; rotation=0., repeatx=0., repeaty=0., shiftx=0.)

Creates a plane of infinitely repeating cylinders in the y-z plane.

# Arguments
- radius: cylinder radius in micrometer
- rotation: angle of cylinders with respect to the z-axis
- repeatx: distance between repeating cylinders along x-axis
- repeaty: distance between repeating cylinders in y-z plane
- shiftx: distance to shift first cylinder from origin along x-axis
"""
function cylinder_plane(radius :: Real; rotation=0., repeatx=0., repeaty=0., shiftx=0.)
    Transformed(
        Repeated(
            Cylinder(radius),
            SA[repeatx, repeaty, 0]
        ),
        CoordinateTransformations.AffineMap(
            Rotations.AngleAxis(Float(deg2rad(rotation)), 1, 0, 0),
            SA[shiftx, 0., 0.]
        )
    )
end
