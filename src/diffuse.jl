struct Movement
    origin :: SVector{3, Real}
    destination :: SVector{3, Real}
    timestep :: Real
end

abstract type Obstruction end

draw_step(diffusivity :: Real, timestep :: Real) = sqrt(timestep * diffusivity) * @SVector randn(3)
draw_step(current_pos :: SVector, diffusivity :: Real, timestep :: Real) = Movement(
    current_pos,
    current_pos + draw_step(diffusivity, timestep),
    timestep
)
function draw_step(current_pos :: SVector, diffusivity :: Real, timestep :: Real, geometry :: Vector{T})  where T <: Obstruction
    correct_collisions(
        draw_step(current_pos, diffusivity, timestep),
        geometry
    )
end

function correct_collisions(to_try :: Movement, geometry :: Vector{T}) where T <: Obstruction
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
    Collision(distance, normal) = new(distance == 0. ? 0. : (distance - eps(distance) * 10.), normal)
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

struct Transformed <: Obstruction
    obstruction :: Obstruction
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