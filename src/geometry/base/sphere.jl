"""
    Sphere([radius[, location]])

Creates a hollow sphere with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
Generate spheres using [`spheres`](@ref).
"""
struct Sphere <: BaseObstruction{3}
    radius :: Float
    id :: UUID
    Sphere(radius) = new(Float(radius), uuid1())
end
Base.copy(s::Sphere) = Sphere(s.radius)

"""
    spheres(radii; positions=[[0, 0, 0]], repeats=[Inf, Inf, Inf], rotation=I(3))

Creates one or more [`Sphere`](@ref)s with given radius (or vector of `radii`).
The `positions`, `repeats`, and `rotation` control the sphere positions and is explained in 
more detail in [Defining the geometry](@ref).
"""
function spheres(args...; kwargs...)
    TransformObstruction(Sphere, args...; kwargs...)
end

isinside(sphere::Sphere, pos::PosVector) = norm(pos) <= sphere.radius
BoundingBox(s::Sphere) = BoundingBox([-s.radius, -s.radius, -s.radius], [s.radius, s.radius, s.radius])

function detect_collision(movement :: Movement, sphere :: Sphere, previous=empty_collision) 
    inside = previous.id != sphere.id ? -1 : previous.index
    sphere_collision(movement, sphere, inside)
end

function sphere_collision(movement :: Movement{N}, obstruction::Obstruction, inside_index::Int) where {N}
    radius = obstruction.radius
    origin = movement.origin
    destination = movement.destination

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
        movement.orientations,
        obstruction.id,
        index=Int(inside)
    )
end

"""
    random_spheres(target_density; repeats, distribution=Distributions.Gamma, mean_radius=1., variance_radius=0.5, max_iter=1000, rotation=I(3))

Generate infinitely repeating box with non-overlapping spheres.

A box with the size of `repeats` will be filled with spheres for a total volume density of `target_density`.
The sphere radii will be drawn from the selected `distribution` (if not set, a Gamma distribution is used with given `mean_radius` and `var_radius`).
An error is raised if no solution for non-overlapping spheres is found.
Other sphere parameters (besides `radii`, `positions`, and `repeats`) are identical as in `mr.spheres`.
"""
function random_spheres(target_density; repeats, distribution=nothing, mean_radius=1., variance_radius=0.5, max_iter=1000, kwargs...)
    (positions, radii) = random_positions_radii(repeats, target_density, 3; distribution=distribution, mean=mean_radius, variance=variance_radius, max_iter=max_iter)
    spheres(radii; positions=positions, repeats=repeats, kwargs...)
end