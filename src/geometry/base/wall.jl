"""
    Wall()

Wall stretching to infinite along two dimensions.
Generate walls using [`walls`](@ref).
"""
struct Wall <: BaseObstruction{1}
    id :: UUID
    Wall() = new(uuid1())
end

Base.copy(w::Wall) = Wall()
isinside(wall::Wall, pos::PosVector) = false
BoundingBox(wall::Wall) = BoundingBox(SA[0], SA[0])

"""
    walls(positions=0, repeats=Inf, rotation=I(3))

Creates one or more [`Wall`](@ref)s.
The `positions`, `repeats`, and `rotation` control the wall position and orientation and is explained in 
more detail in [Defining the geometry](@ref).
"""
walls(;kwargs...) = TransformObstruction(Wall(); kwargs...)

function detect_collision(movement :: Movement{1}, wall :: Wall, previous=empty_collision)
    if previous.id == wall.id
        return empty_collision
    end
    origin = movement.origin[1]
    destination = movement.destination[1]
    if origin * destination > 0
        return empty_collision
    end
    total_length = abs(origin - destination)
    Collision(
        abs(origin) / total_length,
        origin < 0 ? SA[-1, 0, 0] : SA[1, 0, 0],
        movement.orientations,
        wall.id
    )
end

