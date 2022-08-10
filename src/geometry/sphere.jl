"""
    Sphere([radius[, location]])

Creates a hollow sphere with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
"""
struct Sphere <: Obstruction
    radius :: Float
    Sphere(radius) = new(Float(radius))
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

isinside(pos::PosVector, sphere::Sphere) = norm(pos) <= sphere.radius

detect_collision(movement :: Movement, sphere :: Sphere, origin::PosVector) = sphere_collision(movement.origin, movement.destination, sphere.radius, origin)

function sphere_collision(origin :: PosVector, destination :: PosVector, radius :: Float, start::PosVector)
    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    if (
        (origin[1] > radius && destination[1] > radius) ||
        (origin[1] < -radius && destination[1] < -radius) ||
        (origin[2] > radius && destination[2] > radius) ||
        (origin[2] < -radius && destination[2] < -radius) ||
        (origin[3] > radius && destination[3] > radius) ||
        (origin[3] < -radius && destination[3] < -radius)
    )
        return nothing
    end
    inside = norm(start) < radius
    diff = destination - origin
    a = sum(diff .* diff)
    b = sum(2 .* origin .* diff)
    c = sum(origin .* origin)
    determinant = b ^ 2 - 4 * a * (c - radius ^ 2)
    if determinant < 0
        return nothing
    end
    sd = sqrt(determinant)
    ai = inv(a)

    # first try solution with lower s
    solution = (inside ? (-b + sd) : (-b - sd)) * 0.5 * ai
    if solution > 1 || solution <= 0
        return nothing
    end
    point_hit = solution * destination + (1 - solution) * origin
    return Collision(
        solution,
        inside ? -point_hit : point_hit
    )
end
