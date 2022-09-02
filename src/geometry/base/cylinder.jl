"""
    Cylinder(radius; chi_I=-0.1, chi_A=-0.1, g_ratio=1)

Creates a hollow cylinder with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
Generate cylinders using [`cylinders`](@ref).
See [Myelinated cylinders](@ref Myelinated_cylinders) for an explanation of the myelin sheath.
"""
struct Cylinder <: BaseObstruction{2}
    radius :: Float
    id :: UUID
    g_ratio :: Float
    chi_I :: Float
    chi_A :: Float
    internal_field :: Float
    external_field :: Float
    function Cylinder(radius; chi_I=-0.1, chi_A=-0.1, g_ratio=1.)
        if isone(g_ratio)
            internal_field, external_field = zero(Float), zero(Float)
        else
            internal_field = -0.75 * chi_A * log(g_ratio)
            external_field = 2 * (chi_I + chi_A / 4) * (1 - g_ratio^2) / (1 + g_ratio)^2 * radius^2
        end
        new(Float(radius), uuid1(), g_ratio, chi_I, chi_A, internal_field, external_field)
    end
end

Base.copy(c::Cylinder) = Cylinder(c.radius; chi_I=c.chi_I, chi_A=c.chi_A, g_ratio=c.g_ratio)
isinside(cyl::Cylinder, pos::SVector{2, Float}) = (pos[1] * pos[1] + pos[2] * pos[2]) <= (cyl.radius * cyl.radius)
BoundingBox(c::Cylinder) = BoundingBox([-c.radius, -c.radius], [c.radius, c.radius])

function total_susceptibility(c::Cylinder)
    r_outer = 2 * c.radius / (1 + c.g_ratio)
    r_inner = c.g_ratio * r_outer
    chi = c.chi_I + c.chi_A / 4
    2 * π * chi * (r_outer^2 - r_inner^2)
end

"""
    cylinders(radii; g_ratio=1, chi_I=-0.1, chi_A=-0.1, positions=[0, 0], repeats=[Inf, Inf], rotation=I(3)

Creates one or more [`Cylinder`](@ref)s with given radius (or vector of `radii`).
[Myelinated cylinders](@ref Myelinated_cylinders) can be created by setting the `g_ratio` to a different value that 1.
All parameters can be either a single value or a vector of values.
The `positions`, `repeats`, and `rotation` control the cylinder position and orientation and is explained in 
more detail in [Defining the geometry](@ref).
"""
function cylinders(args...; kwargs...)
    TransformObstruction(Cylinder, args...; kwargs...)
end

function detect_collision(movement :: Movement{2}, cylinder :: Cylinder, previous=empty_collision)
    inside = previous.id != cylinder.id ? -1 : previous.index
    sphere_collision(movement.origin, movement.destination, cylinder, inside)
end


"""
    off_resonance(cylinder, position, b0_field)

Computed by the hollow cylinder fiber model from [Wharton_2012](@cite).
"""
function off_resonance(cylinder::Cylinder, position::SVector{2, Float}, b0_field::SVector{2, Float})
    if iszero(cylinder.internal_field) && iszero(cylinder.external_field)
        return zero(Float)
    end
    rsq = position[1] * position[1] + position[2] * position[2]
    sin_theta_sq = 1 - b0_field[1] * b0_field[1] - b0_field[2] + b0_field[2]
    if rsq < cylinder.radius^2
        return cylinder.internal_field * sin_theta_sq
    else
        cos2 = (b0_field[1] * position[1] + b0_field[2] * position[2])^2 / rsq
        # cos 2 phi = cos^2 phi - sin^2 phi = 1 - 2 cos^2 phi
        return cylinder.external_field * sin_theta_sq * (2 * cos2 - 1) / rsq
    end
end

function lorentz_off_resonance(cylinder::Cylinder, position::SVector{2, Float}, b0_field::SVector{2, Float}, repeat_dist::SVector{2, Float}, radius::Float, nrepeats::SVector{2, Int})
    field = zero(Float)
    if iszero(cylinder.internal_field) && iszero(cylinder.external_field)
        return field
    end
    lorentz_radius_sq = radius * radius
    sin_theta_sq = b0_field[1] * b0_field[1] + b0_field[2] + b0_field[2]
    if iszero(sin_theta_sq)
        return field
    end
    for i in -nrepeats[1]:nrepeats[1]
        xshift = i * repeat_dist[1]
        p1 = position[1] + xshift
        for j in -nrepeats[2]:nrepeats[2]
            yshift = j * repeat_dist[2]
            p2 = position[2] + yshift
            rsq = p1 * p1 + p2 * p2
            if rsq > lorentz_radius_sq
                continue
            end
            if i == 0 && j == 0 && rsq < cylinder.radius^2
                field += cylinder.internal_field * sin_theta_sq
            else
                cos2 = (b0_field[1] * p1 + b0_field[2] * p2)^2 / rsq
                field += cylinder.external_field * sin_theta_sq * (2 * cos2 - 1) / rsq
            end
        end
    end
    return field
end

"""
    random_cylinders(target_density; repeats, distribution=Distributions.Gamma, mean_radius=1., variance_radius=0.5, max_iter=1000, g_ratio=1., chi_I=-0.1, chi_A=-0.1, rotation=I(3))

Generate infinitely repeating box with non-overlapping cylinders.

A rectangle with the size of `repeats` will be filled with cylinders for a total surface density of `target_density`.
The cylinder radii will be drawn from the selected `distribution` (if not set, a Gamma distribution is used with given `mean_radius` and `var_radius`).
An error is raised if no solution for non-overlapping cylinders is found.
Other cylinder parameters (besides `radii`, `shifts`, and `repeats`) are identical as in `mr.cylinders`.
"""
function random_cylinders(target_density; repeats, distribution=nothing, mean_radius=1., variance_radius=0.5, max_iter=1000, kwargs...)
    (positions, radii) = random_positions_radii(repeats, target_density, 2; distribution=distribution, mean=mean_radius, variance=variance_radius, max_iter=max_iter)
    cylinders(radii; positions=positions, repeats=repeats, kwargs...)
end