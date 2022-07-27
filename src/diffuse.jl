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

struct Wall <: Obstruction
    normal :: SVector{3, Real}
    offset :: Real
    function Wall(normal, offset)
        n = norm(normal)
        new(normal / n, offset * n)
    end
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
    origin = movement.origin ⋅ wall.normal - wall.offset
    destination = movement.destination ⋅ wall.normal - wall.offset
    if origin * destination >= 0
        return nothing
    end
    total_length = abs(origin - destination)
    Collision(
        abs(origin) / total_length,
        wall.normal
    )
end

struct Sphere <: Obstruction
    radius :: Real
    location :: SVector{3, Real}
end

function detect_collision(movement :: Movement, sphere :: Sphere)
    origin = movement.origin .- sphere.location
    destination = movement.destination .- sphere.location

    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    a = sum((destination .- origin) .^ 2)
    b = sum(2 .* origin .* (destination .- origin))
    c = sum(origin .* origin)
    determinant = b ^ 2 - 4 * a * (c - sphere.radius ^ 2)
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
    orientation :: SVector{3, Real}
    location :: SVector{3, Real}
    function Cylinder(radius :: Real, orientation :: SVector, location :: SVector)
        n_orientation = orientation / norm(orientation)
        rel_offset = location - (location ⋅ n_orientation) * n_orientation
        new(radius, n_orientation, rel_offset)
    end
end

function Cylinder(radius :: Real, sym :: Symbol, offset :: SVector)
    orientation = Dict(
        :x => SA_F64[1., 0., 0.],
        :y => SA_F64[0., 1., 0.],
        :z => SA_F64[0., 0., 1.],
    )
    Cylinder(radius, orientation[sym], offset)
end

function normed_offset(position :: SVector, cylinder :: Cylinder)
    offset = position .- cylinder.location
    offset .- (offset ⋅ cylinder.orientation) .* cylinder.orientation
end

function detect_collision(movement :: Movement, cylinder :: Cylinder)
    sphere = Sphere(cylinder.radius, cylinder.location)
    rel_movement = Movement(
        normed_offset(movement.origin, cylinder),
        normed_offset(movement.destination, cylinder),
        movement.timestep
    )
    detect_collision(rel_movement, sphere)
end