"""
    Wall([[normal,] offset])

Infinitely large wall with an orientation given by `normal` (default in the x-direction).
The normal can also be defined using :x, :y, or :z, to point in that cardinal direction.
The offset of the wall from the origin along this `normal` is given by `offset` (so that `offset .* normal` is on the wall).
"""
struct Wall <: Obstruction
    id :: UUID
    Wall() = new(uuid1())
end

Wall(offset :: Real) = iszero(offset) ? Wall() : Transformed(Wall(), CoordinateTransformations.Translation(Float(offset), zero(Float), zero(Float)))

isinside(pos::PosVector, wall::Wall) = false
BoundingBox(wall::Wall) = BoundingBox(SA[0, -Inf, -Inf], SA[0, Inf, Inf])


function Wall(normal :: AbstractVector, offset :: Real)
    normal = PosVector(normal)
    offset = Float(offset)
    n = normal ./ norm(normal)
    shifted = Wall(offset * norm(normal))
    if isapprox(n, SA[1, 0, 0], atol=1e-10)
        return shifted
    end
    rot_axis = cross(SA[1, 0, 0], n)
    rot_angle = acos(n[1])
    Transformed(shifted, CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Wall(sym :: Symbol, offset :: Real)
    direction = Dict(
        :x => SA[1., 0., 0.],
        :y => SA[0., 1., 0.],
        :z => SA[0., 0., 1.],
    )
    Wall(direction[sym], offset)
end

function detect_collision(movement :: Movement, wall :: Wall, previous=empty_collision)
    if previous.id == wall.id
        return empty_collision
    end
    origin = movement.origin[1]
    destination = movement.destination[1]
    if origin * destination >= 0
        return empty_collision
    end
    total_length = abs(origin - destination)
    Collision(
        abs(origin) / total_length,
        origin < 0 ? SA[-1, 0, 0] : SA[1, 0, 0],
        wall.id
    )
end

