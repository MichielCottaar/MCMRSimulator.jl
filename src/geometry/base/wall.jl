"""
    Wall()

Wall stretching to infinite along two dimensions.
Generate walls using [`walls`](@ref).
"""
struct Wall <: BaseObstruction{1}
    properties :: ObstructionProperties
end
Wall(; kwargs...) = Wall(ObstructionProperties(; kwargs...))

isinside(wall::Wall, pos::PosVector) = false
BoundingBox(wall::Wall) = BoundingBox(SA[0], SA[0])

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
    if (id(previous) == id(wall)) && (abs(origin) < 1e-12)
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

