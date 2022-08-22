"""
    Cylinder([radius[, orientation[, location]; chi_I=-0.1, chi_A=-0.1, g_ratio=1]])

Creates a hollow cylinder with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
The orientation of the cylinder (default: z-direction) can be given as a symbol of the cardinal orientation (:x, :y, or :z) or as a length-3 vector.

# Myelin sheath
The cylinder can be wrapped in a myelin sheath by setting `g_ratio` to a different value from 1.
The myelin sheath will be infinitely thin, when considering collisions,
however it will generate an off-resonance field determined by the myelin's isotropic magnetic susceptibility (`chi_I`),
the anisotropic magnetic susceptibility (`chi_A`), and the g-ratio (`g_ratio`).
Both `chi_I` and `chi_A` are given in ppm and set to the value from Wharton & Bowtell (2012).
"""
struct Cylinder <: Obstruction
    radius :: Float
    id :: UUID
    internal_field :: Float
    external_field :: Float
    function Cylinder(radius; chi_I=-0.1, chi_A=-0.1, g_ratio=1.)
        if isone(g_ratio)
            internal_field, external_field = zero(Float), zero(Float)
        else
            internal_field = -0.75 * chi_A * log(g_ratio)
            external_field = (chi_I + chi_A / 4) * (1 - g_ratio^2) * (1 + g_ratio^2) / 4 * radius^2
        end
        new(Float(radius), uuid1(), internal_field, external_field)
    end
end

isinside(pos::PosVector, cyl::Cylinder) = (pos[1] * pos[1] + pos[2] * pos[2]) <= (cyl.radius * cyl.radius)
BoundingBox(c::Cylinder) = BoundingBox([-c.radius, -c.radius, -Inf], [c.radius, c.radius, Inf])

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real}; kwargs...)
    radius = Float(radius)
    o = SVector{3, Float}(orientation / norm(orientation))
    if isapprox(o, SA[0, 0, 1], atol=1e-10)
        return Cylinder(radius; kwargs...)
    end
    rot_axis = cross(SA[0, 0, 1], o)
    rot_angle = acos(o[end])
    Transformed(Cylinder(radius; kwargs...), CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real}, location :: AbstractVector{<:Real}; kwargs...)
    c = Cylinder(radius, orientation; kwargs...)
    location = SVector{3, Float}(location)
    return all(iszero.(location)) ? c : Transformed(c, CoordinateTransformations.Translation(location...))
end

function Cylinder(radius :: Real, sym :: Symbol, args...; kwargs...)
    orientation = Dict(
        :x => SA[1., 0., 0.],
        :y => SA[0., 1., 0.],
        :z => SA[0., 0., 1.],
    )
    Cylinder(radius, orientation[sym], args...; kwargs...)
end

Cylinder(;radius=1., orientation=[0., 0, 1], position=[0., 0, 0]) = Cylinder(radius, orientation, position)

function detect_collision(movement :: Movement, cylinder :: Cylinder, previous=empty_collision)
    select(a) = SA[a[1], a[2]]
    inside = previous.id != cylinder.id ? -1 : previous.index
    sphere_collision(select(movement.origin), select(movement.destination), cylinder, inside)
end


"""
    cylinder_plane(radius; rotation=0., repeatx=0., repeaty=0., shiftx=0.)

Creates a plane of infinitely repeating cylinders in the y-z plane.

# Arguments
- radius: cylinder radius in micrometer
- rotation: angle of cylinders with respect to the z-axis
- repeatx: distance between repeating cylinders along x-axis
- repeaty: distance between repeating cylinders in y-z plane
- shiftx: distance to shift first cylinder from origin along x-axis
"""
function cylinder_plane(radius :: Real; rotation=0., repeatx=0., repeaty=0., shiftx=0.)
    Transformed(
        Repeated(
            Cylinder(radius),
            SA[repeatx, repeaty, 0]
        ),
        CoordinateTransformations.AffineMap(
            Rotations.AngleAxis(Float(deg2rad(rotation)), 1, 0, 0),
            SA[shiftx, 0., 0.]
        )
    )
end


"""
    off_resonance(position, cylinder, b0_field)

Computed by the hollow cylinder fiber model from Wharton & Bowtell (2012).
"""
function off_resonance(cylinder::Cylinder, position::PosVector, b0_field::PosVector)
    if iszero(cylinder.internal_field) && iszero(cylinder.external_field)
        return zero(Float)
    end
    println("computing field")
    println(position)
    println(b0_field)
    rsq = position[1] * position[1] + position[2] * position[2]
    cos_theta_sq = b0_field[3]^2  # theta is the angle between the b0_field and the cylinder orientation
    sin_theta_sq = 1 - cos_theta_sq
    if rsq < cylinder.radius^2
        return cylinder.internal_field * sin_theta_sq
    else
        cos2 = (b0_field[1] * position[1] + b0_field[2] * position[2])^2 / rsq
        # cos 2 phi = cos^2 phi - sin^2 phi = 1 - 2 cos^2 phi
        return cylinder.external_field * sin_theta_sq * (2 * cos2 - 1) / rsq
    end
end