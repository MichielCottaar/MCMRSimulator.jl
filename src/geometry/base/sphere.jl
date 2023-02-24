"""
    Sphere([radius[, location]])

Creates a hollow sphere with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
Generate spheres using [`spheres`](@ref).
"""
struct Sphere <: BaseObstruction{3}
    radius :: Float
    properties :: ObstructionProperties
end
Sphere(radius; kwargs...) = Sphere(Float(radius), ObstructionProperties(; kwargs...))

"""
    spheres(radii; positions=[[0, 0, 0]], repeats=[Inf, Inf, Inf], rotation=I(3))

Creates one or more [`Sphere`](@ref)s with given radius (or vector of `radii`).
The `positions`, `repeats`, and `rotation` control the sphere positions and is explained in 
more detail in [Defining the geometry](@ref).
Additional keyword arguments are available to set generic obstruction settings as described in [`ObstructionProperties`](@ref).
"""
function spheres(args...; kwargs...)
    TransformObstruction(Sphere, args...; kwargs...)
end

isinside(sphere::Sphere, pos::PosVector, stuck_to::Collision) = Int(collided(sphere, stuck_to) ? stuck_to.inside : norm(pos) <= sphere.radius)
BoundingBox(s::Sphere) = BoundingBox{3}(s.radius)

function detect_collision(movement :: Movement, sphere :: Sphere, previous=empty_collision) 
    inside = !collided(sphere, previous) ? -1 : Int(previous.inside)
    sphere_collision(movement, sphere, inside)
end

function sphere_collision(movement :: Movement{N}, obstruction::Obstruction, inside_index::Int) where {N}
    radius = obstruction.radius
    origin = movement.origin
    destination = movement.destination

    rsq_origin = sum(origin .* origin)
    inside = inside_index == -1 ? rsq_origin < radius * radius : Bool(inside_index)
    diff = destination - origin

    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    a = sum(diff .* diff)
    b = sum(2 .* origin .* diff)
    c = rsq_origin
    determinant = b * b - 4 * a * (c - radius ^ 2)
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
        normal = SA[point_hit[1], point_hit[2], 0.]
    else
        normal = point_hit
    end
    return Collision(
        solution,
        inside ? -normal : normal,
        ObstructionProperties(obstruction),
        inside=inside
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

size_scale(sphere::Sphere) = sphere.radius

"""
    random_on_sphere()

Draws a random orientation on the unit sphere as a length-3 vector

The z-orientation is drawn a random number between -1 and 1.
The angle in the x-y plane is drawn as a random number between 0 and 2π.
This results in an unbiased random distribution across the sphere.
"""
function random_on_sphere()
    z = rand(Float) * 2 - 1
    r = sqrt(1 - z*z)
    theta = rand(Float) * Float(2 * π)
    (s, c) = sincos(theta)
    return SA[
        r * s,
        r * c,
        z
    ]
end

function random_surface_positions(sphere::Sphere, total_density::Number)
    nspins = Int(floor(total_density * 4π * sphere.radius^2 + rand()))
    normals = [random_on_sphere() for _ in 1:nspins]
    positions = -normals .* sphere.radius
    return (positions, normals, zeros(Int, nspins))
end