"""
Defines a finite plane in the 3D space used for plotting.

# Constructor
    PlotPlane(normal::PosVector, position::PosVector; sizex=Inf, sizey=Inf, ngrid=100)

Arguments:
- `normal`: length-3 vector with the orientation perpendicular to the plane (default to z-direction).
- `position`: position of plane as a length-3 vector (defaults to origin).
- `sizex`: size of the plane in the x-direction (before rotating to `normal`).
- `sizey`: size of the plane in the y-direction (before rotating to `normal`).
- `ngrid`: number of grid elements to split the plane up into for plotting.

# Spin projection onto plane
See [`project`](@ref) for details on how spins are projected onto the `PlotPlane`.
"""
struct PlotPlane
    transformation :: CoordinateTransformations.Transformation
    sizex :: Real
    sizey :: Real
    ngrid :: Int
end


function PlotPlane(
    normal :: AbstractVector{<:Real}=SA[0, 0, 1], 
    position :: AbstractVector{<:Real}=SA[0, 0, 0];
    sizex::Real=Inf, sizey::Real=Inf,
    ngrid=100,
)
    normal = SVector{3}(normal)
    position = SVector{3}(position)
    if normal ≈ SA[0, 0, 1]
        transform = CoordinateTransformations.Translation(position)
    else
        rot_axis = cross(SA[0, 0, 1], normal)
        rot_angle = acos(normal[3] / norm(normal))
        transform = CoordinateTransformations.AffineMap(
            Rotations.AngleAxis(rot_angle, rot_axis...),
            position
        )
    end
    PlotPlane(CoordinateTransformations.inv(transform), sizex, sizey, ngrid)
end

"""
    project(plot_plane, position)
    project(plot_plane, snapshot)

Transforms the `position` (length-3 vector) or [`Snapshot`](@ref) to a space, where the [`PlotPlane`](@ref) lies in the x-y-plane centered on origin.
"""
function project(pp::PlotPlane, pos::PosVector)
    base = pp.project(pos)
    sizex, sizey = pp.sizex, pp.sizey
    correct = [isfinite(sizex) ? sizex : 0., isfinite(sizey) ? sizey : 0., 0.] / 2.
    mod.(base .+ correct, (sizex, sizey, Inf)) .- correct
end

project(pp::PlotPlane, spin::Spin) = project(pp, mr.position(spin))

function project(pp::PlotPlane, snap::Snapshot)
    Snapshot(
        [Spin(project(pp, s), s.orientation) for s in snap.spins],
        snap.time
    )
end


"""
    project_on_grid(plot_plane, snap)

Spins from the [`Snapshot`](@ref) are projected onto the grid defined by [`PlotPlane`](@ref) in two ways:
- along the normal spins are projected onto the plane from infinitely far (TODO: give finite extent)
- in the other directions any spins are projected onto the plane using mod(position[1], `sizex`) and mod(position[2], `sizey`).
    This assumes that the geometry and field repeats itself ad infinitum beyond the `PlotPlane` (TODO: allow this assumption to be turned off).
In effect, this means that all spins are projected onto the `PlotPlane`. The average spin orientation in each grid cell is returned.
"""
function project_on_grid(pp::PlotPlane, snap::Snapshot)
    on_plane = transform(pp, snap)
    positions = [s.position[1:2] for s in on_plane]
    xrange = extrema(v->v[1], positions)
    yrange = extrema(v->v[2], positions)

    res = zeros(MVector{3}, pp.ngrid, pp.ngrid)
    hits = zeros(Int, pp.ngrid, pp.ngrid)
    for spin in on_plane
        relx = (spin.position[1] - xrange[1])/(xrange[2] - xrange[1])
        rely = (spin.position[2] - yrange[1])/(yrange[2] - yrange[1])
        x_index = min.(Int(floor(relx * (pp.ngrid - 1))), pp.ngrid - 2) + 1
        y_index = min.(Int(floor(rely * (pp.ngrid - 1))), pp.ngrid - 2) + 1
        vs = vector(spin.orientation)
        for (xi, yi) in [
            (x_index, y_index),
            (x_index+1, y_index),
            (x_index, y_index+1),
            (x_index+1, y_index+1),
        ]
            res[xi, yi] += vs
            hits[xi, yi] += 1
        end
    end
    (range(xrange..., pp.ngrid), range(yrange..., pp.ngrid), vector2spin.(res ./ hits))
end

