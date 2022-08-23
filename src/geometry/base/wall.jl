"""
    Wall([[normal,] offset])

Infinitely large wall with an orientation given by `normal` (default in the x-direction).
The normal can also be defined using :x, :y, or :z, to point in that cardinal direction.
The offset of the wall from the origin along this `normal` is given by `offset` (so that `offset .* normal` is on the wall).
"""
struct Wall <: BaseObstruction{1}
    id :: UUID
    Wall() = new(uuid1())
end

Base.copy(w::Wall) = Wall()
isinside(pos::PosVector, wall::Wall) = false
BoundingBox(wall::Wall) = BoundingBox(SA[0], SA[0])

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
        wall.id
    )
end

