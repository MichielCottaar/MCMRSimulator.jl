"""
    Cylinder([radius[, orientation[, location]]])

Creates a hollow cylinder with a radius of `radius` micrometer (default 1 micrometer) at the given `location` (default: origin).
The orientation of the cylinder (default: z-direction) can be given as a symbol of the cardinal orientation (:x, :y, or :z) or as a length-3 vector.
"""
struct Cylinder <: Obstruction
    radius :: Float
    id :: UUID
    Cylinder(radius) = new(Float(radius), uuid1())
end

isinside(pos::PosVector, cyl::Cylinder) = (pos[1] * pos[1] + pos[2] * pos[2]) <= (cyl.radius * cyl.radius)
BoundingBox(c::Cylinder) = BoundingBox([-c.radius, -c.radius, -Inf], [c.radius, c.radius, Inf])

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real})
    radius = Float(radius)
    o = SVector{3, Float}(orientation / norm(orientation))
    if isapprox(o, SA[0, 0, 1], atol=1e-10)
        return Cylinder(radius)
    end
    rot_axis = cross(SA[0, 0, 1], o)
    rot_angle = acos(o[end])
    Transformed(Cylinder(radius), CoordinateTransformations.LinearMap(Rotations.AngleAxis(rot_angle, rot_axis...)))
end

function Cylinder(radius :: Real, orientation :: AbstractVector{<:Real}, location :: AbstractVector{<:Real})
    c = Cylinder(radius, orientation)
    location = SVector{3, Float}(location)
    return all(iszero.(location)) ? c : Transformed(c, CoordinateTransformations.Translation(location...))
end

function Cylinder(radius :: Real, sym :: Symbol, offset :: AbstractVector{<:Real})
    orientation = Dict(
        :x => SA[1., 0., 0.],
        :y => SA[0., 1., 0.],
        :z => SA[0., 0., 1.],
    )
    Cylinder(radius, orientation[sym], offset)
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
