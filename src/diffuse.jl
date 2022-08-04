"""
Trajectory of a single spin through space.

The trajectory can be accessed through indexing at a specific time (`trajectory[time]`) 
or through iterations (`for (time, position) in trajectory`).
The trajectory will be automatically extended as needed for ever.
"""
struct StepTrajectory{T <: AbstractFloat} 
    positions :: Vector{SVector{3, T}}
    timestep :: Real
    diffusivity :: Field
    geometry :: Obstructions
    function StepTrajectory(
        origin::SVector{3, T}, timestep::Real, diffusivity :: Field, geometry :: Obstructions
        ) where T <: AbstractFloat
        new{T}([origin], timestep, diffusivity, geometry)
    end
end

Base.iterate(iter::StepTrajectory) = ((0., iter.positions[1]), 1)
Base.iterate(iter::StepTrajectory, state::Int) = ((iter.timestep * state, iter[state + 1]), state + 1)
Base.eltype(::Type{StepTrajectory{T}}) where {T} = Tuple{Int, SVector{3, T}}
Base.IteratorSize(::Type{StepTrajectory}) = Base.IsInfinite()
Base.firstindex(st :: StepTrajectory) = 1

function Base.getindex(st::StepTrajectory, time::Real)
    if isa(st.diffusivity, ZeroField)
        return st.positions[1]
    end
    index = Int(div(time, st.timestep, RoundNearest)) + 1
    while index > length(st.positions)
        push!(st.positions, draw_step(
            st.positions[end], st.diffusivity(st.positions[end]), 
            st.timestep, st.geometry
            ))
    end
    st.positions[index]
end


function random_on_sphere()
    z = rand() * 2. - 1.
    r = sqrt(1. - z*z)
    theta = rand() * 2 * π
    return SA_F64[
        r * sin(theta),
        r * cos(theta),
        z
    ]
end

draw_step(diffusivity :: Real, timestep :: Real) = sqrt(timestep * diffusivity) * @SVector randn(3)
draw_step(current_pos :: SVector, diffusivity :: Real, timestep :: Real) = current_pos .+ draw_step(diffusivity, timestep)
function draw_step(current_pos :: SVector, diffusivity :: Real, timestep :: Real, geometry :: Obstructions)
    new_pos = draw_step(current_pos, diffusivity, timestep)
    if geometry == Obstruction[]
        return new_pos
    end
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


struct Collision
    distance :: Real
    normal :: SVector{3, Real}
    Collision(distance, normal) = new(distance * (1. - eps(typeof(distance))^0.75), normal)
end


function detect_collision(movement :: Movement, obstructions :: Vector{T}) where T <: Obstruction
    collision = nothing
    for o in obstructions
        c_new = detect_collision(movement, o)
        if !isnothing(c_new) && (isnothing(collision) || c_new.distance < collision.distance)
            collision = c_new
        end
    end
    collision
end



struct Repeated <: Obstruction
    obstruction :: Obstructions
    repeats :: SVector{3, Real}
    function Repeated(obstruction, repeats)
        rs = SVector([iszero(r) ? Inf : abs(r) for r in repeats])
        any(isfinite.(rs)) ? new(obstruction, rs) : obstruction
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
            repeat.obstruction
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

"Computes all interactions with a 1x1x1 grid"
ray_grid_intersections(origin :: SVector, destination :: SVector) = Channel() do c
    direction = destination .- origin
    within_voxel = mod.(origin, 1)
    all_next_hits = MVector{3, Float64}([(d > 0 ? 1 - w : w) / abs(d) for (d, w) in zip(direction, within_voxel)])
    prev_pos = origin
    prev_time = 0.
    current_voxel = MVector{3, Int}(Int.(floor.(origin)))
    while true
        dimension = argmin(all_next_hits)
        time_to_hit = all_next_hits[dimension]
        next_time = prev_time + time_to_hit
        if next_time > 1.
            push!(c, (current_voxel, prev_time, prev_pos .- current_voxel, 1., destination .- current_voxel))
            break
        end
        next_pos = origin .+ direction .* next_time
        for dim in 1:3
            if dim == dimension
                push!(c, (SVector(current_voxel), prev_time, prev_pos .- current_voxel, next_time, next_pos .- current_voxel))
                prev_time = next_time
                prev_pos = next_pos
                all_next_hits[dim] = 1 / abs(direction[dim])
                current_voxel[dim] += sign(direction[dim])
            else
                all_next_hits[dim] -= time_to_hit
            end
        end
    end
end

struct Transformed <: Obstruction
    obstruction :: Obstructions
    transform :: CoordinateTransformations.Transformation
    inverse :: CoordinateTransformations.Transformation
    function Transformed(obstruction :: Obstruction, transform :: CoordinateTransformations.Transformation)
        if isa(obstruction, Transformed)
            return Transformed(obstruction.obstruction, transform ∘ obstruction.transform )
        else
            return new(obstruction, transform, inv(transform))
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
        transform.obstruction
    )
    if isnothing(c)
        return c
    end
    Collision(
        c.distance,
        transform.transform(c.normal) .- transform.transform(zero(SVector{3, Real}))
    )
end

struct Wall <: Obstruction
end

Wall(offset :: Real) = offset == 0 ? Wall() : Transformed(Wall(), CoordinateTransformations.Translation(offset, 0., 0.))


function Wall(normal :: SVector, offset :: Real)
    n = normal ./ norm(normal)
    shifted = Wall(offset * norm(normal))
    if isapprox(n, SA_F64[1, 0, 0], atol=1e-10)
        return shifted
    end
    rot_axis = cross(SA_F64[1, 0, 0], n)
    rot_angle = acos(n[1])
    Transformed(shifted, CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Wall(sym :: Symbol, offset :: Real)
    direction = Dict(
        :x => SA_F64[1., 0., 0.],
        :y => SA_F64[0., 1., 0.],
        :z => SA_F64[0., 0., 1.],
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
        SA_F64[1, 0, 0]
    )
end

struct Sphere <: Obstruction
    radius :: Real
end

function Sphere(radius :: Real, location :: SVector)
    if all(iszero.(location))
        return Sphere(radius)
    else
        return Transformed(Sphere(radius), CoordinateTransformations.Translation(location))
    end
end

detect_collision(movement :: Movement, sphere :: Sphere) = sphere_collision(movement.origin, movement.destination, sphere.radius)

function sphere_collision(origin :: AbstractVector, destination :: AbstractVector, radius :: Real)
    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    a = sum((destination .- origin) .^ 2)
    b = sum(2 .* origin .* (destination .- origin))
    c = sum(origin .* origin)
    determinant = b ^ 2 - 4 * a * (c - radius ^ 2)
    if determinant < 0
        return nothing
    end

    # first try solution with lower s
    solution = (-b - sqrt(determinant)) / (2 * a)
    if solution > 1
        return nothing
    elseif solution <= 0
        solution = (-b + sqrt(determinant)) / (2 * a)
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

struct Cylinder <: Obstruction
    radius :: Real
end

function Cylinder(radius :: Real, orientation :: SVector)
    o = orientation / norm(orientation)
    if isapprox(o, SA_F64[0, 0, 1], atol=1e-10)
        return Cylinder(radius)
    end
    rot_axis = cross(SA_F64[0, 0, 1], o)
    rot_angle = acos(o[end])
    Transformed(Cylinder(radius), CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Cylinder(radius :: Real, orientation :: SVector, location :: SVector)
    c = Cylinder(radius, orientation)
    return all(iszero.(location)) ? c : Transformed(c, CoordinateTransformations.Translation(location...))
end

function Cylinder(radius :: Real, sym :: Symbol, offset :: SVector)
    orientation = Dict(
        :x => SA_F64[1., 0., 0.],
        :y => SA_F64[0., 1., 0.],
        :z => SA_F64[0., 0., 1.],
    )
    Cylinder(radius, orientation[sym], offset)
end

function detect_collision(movement :: Movement, cylinder :: Cylinder)
    o = SA_F64[movement.origin[1], movement.origin[2], 0.]
    d = SA_F64[movement.destination[1], movement.destination[2], 0.]
    sphere_collision(o, d, cylinder.radius)
end


function cylinder_plane(radius :: Real; rotation=0., repeatx=0., repeaty=0., shiftx=0.)
    Transformed(
        Repeated(
            Cylinder(radius),
            SA_F64[repeatx, repeaty, 0]
        ),
        CoordinateTransformations.AffineMap(
            Rotations.AngleAxis(deg2rad(rotation), 1, 0, 0),
            [shiftx, 0, 0]
        )
    )
end
