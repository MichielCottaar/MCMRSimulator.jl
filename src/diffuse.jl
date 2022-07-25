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
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal .+ direction
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
    Collision(distance, normal) = new(distance == 0. ? 0. : prevfloat(distance), normal)
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