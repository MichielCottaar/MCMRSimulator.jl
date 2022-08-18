"""
    Sphere([radius[, location]])

Creates a hollow sphere with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
"""
struct Sphere <: Obstruction
    radius :: Float
    id :: UUID
    Sphere(radius) = new(Float(radius), uuid1())
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
Sphere(;radius=1., position=[0, 0, 0]) = Sphere(radius, position)

isinside(pos::PosVector, sphere::Sphere) = norm(pos) <= sphere.radius
BoundingBox(s::Sphere) = BoundingBox([-s.radius, -s.radius, -s.radius], [s.radius, s.radius, s.radius])

function detect_collision(movement :: Movement, sphere :: Sphere, previous=empty_collision) 
    inside = previous.id != sphere.id ? -1 : previous.index
    sphere_collision(movement.origin, movement.destination, sphere, inside)
end

function sphere_collision(origin :: SVector{N, Float}, destination :: SVector{N, Float}, obstruction::Obstruction, inside_index::Int) where {N}
    radius = obstruction.radius
    for dim in 1:N
        if (
            (origin[dim] > radius && destination[dim] > radius) ||
            (origin[dim] < -radius && destination[dim] < -radius)
        )
            return empty_collision
        end
    end
    inside = inside_index == -1 ? norm(origin) < radius : Bool(inside_index)
    diff = destination - origin

    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    a = sum(diff .* diff)
    b = sum(2 .* origin .* diff)
    c = sum(origin .* origin)
    determinant = b ^ 2 - 4 * a * (c - radius ^ 2)
    if determinant < 0
        return empty_collision
    end
    sd = sqrt(determinant)
    ai = inv(a)

    # first try solution with lower s
    solution = (inside ? (-b + sd) : (-b - sd)) * 1//2 * ai
    if solution > 1 || solution <= 0
        return empty_collision
    end
    point_hit = solution * destination + (1 - solution) * origin
    if N == 2
        point_hit = SA[point_hit[1], point_hit[2], 0.]
    end
    return Collision(
        solution,
        inside ? -point_hit : point_hit,
        obstruction.id,
        index=Int(inside)
    )
end
