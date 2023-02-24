"""
    Wall()

Wall stretching to infinite along two dimensions.
Generate walls using [`walls`](@ref).
"""
struct Wall <: BaseObstruction{1}
    properties :: ObstructionProperties
end
Wall(; kwargs...) = Wall(ObstructionProperties(; kwargs...))

isinside(wall::Wall, pos::SVector{1, Float}, stuck_to::Collision) = 0
BoundingBox(wall::Wall) = BoundingBox{1}(0)

"""
    walls(positions=0, repeats=Inf, rotation=I(3))

Creates one or more [`Wall`](@ref)s.
The `positions`, `repeats`, and `rotation` control the wall position and orientation and is explained in 
more detail in [Defining the geometry](@ref).
Additional keyword arguments are available to set generic obstruction settings as described in [`ObstructionProperties`](@ref).
"""
walls(;kwargs...) = TransformObstruction(Wall; kwargs...)

function detect_collision(movement :: Movement{1}, wall :: Wall, previous=empty_collision)
    origin = movement.origin[1]
    if collided(wall, previous) && (abs(origin) < 1e-12)
        return empty_collision
    end
    destination = movement.destination[1]
    if origin * destination > 0
        return empty_collision
    end
    total_length = abs(origin - destination)
    Collision(
        abs(origin) / total_length,
        origin < 0 ? SA[-1, 0, 0] : SA[1, 0, 0],
        wall.properties,
        inside=origin > 0
    )
end

size_scale(wall::Wall) = Inf

function random_surface_positions(wall::Wall, total_density::Number)
    nspins = Int(floor(total_density + rand()))
    return (zeros(SVector{1, Float}, nspins), ones(SVector{1, Float}), zeros(Int, nspins))
end
