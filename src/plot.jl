import Makie: Makie, @lift

"""
    color(orient::SpinOrientation; saturation=1.)

Returns a color representing the spin orientation in the transverse (x-y) plane.
Brighter colors have a larger transverse component, so that spins with no transverse component are black.
The actual color encodes the spin orientation.
"""
color(orient::Union{Spin, SpinOrientation}; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))

@Makie.recipe(SequencePlot, seq) do scene
    Makie.Theme(
    )
end

"""
    plot(sequence)
    plot!(sequence)

Creates a visual representation of a [`Sequence`](@ref) diagram.
"""
function Makie.plot!(sp::SequencePlot)
    seq = sp[1]
    times = @lift [p.time for p in $seq.pulses if isa(p, RFPulse)]
    flip_angle = @lift [rad2deg(p.flip_angle) for p in $seq.pulses if isa(p, RFPulse)]
    as_zero = @lift $flip_angle * 0
    max_angle = @lift maximum($flip_angle)

    yval = @lift $flip_angle ./ $max_angle
    Makie.arrows!(sp, times, as_zero, as_zero, yval)

    text = @lift [string(Int(round(a))) for a in $flip_angle]
    text_positions = @lift [(t, 0.05 + y) for (t, y) in zip($times, $yval)]
    Makie.text!(sp, text, position=text_positions, align=(:center, :center))
    #vlines!(sp, TR, color="black")
    sp
end

Makie.plottype(::Sequence) = SequencePlot


@Makie.recipe(SnapshotPlot, snap) do scene
    Makie.Theme(
    )
end


"""
    plot(snapshot)
    plot!(snapshot)

Plots the spin positions in the [`Snapshot`](@ref) in 3D color coded by the spin's orientation (see [`color`](@ref)).
"""
function Makie.plot!(sp::SnapshotPlot)
    snap = sp[1]
    colors = @lift color.($snap.spins)
    pos = @lift [Makie.Point3f(s.position) for s in $snap.spins]
    Makie.meshscatter!(sp, pos, color=colors)
    sp
end

Makie.plottype(::Snapshot) = SnapshotPlot


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
See [`project_on_grid`](@ref) for details on how spins are projected onto the `PlotPlane`.
"""
struct PlotPlane
    transformation :: CoordinateTransformations.Transformation
    sizex :: Real
    sizey :: Real
    ngrid :: Int
end


function PlotPlane(
    normal :: AbstractVector{<:Real}=SA_F64[0, 0, 1], 
    position :: AbstractVector{<:Real}=SA_F64[0, 0, 0];
    sizex::Real=Inf, sizey::Real=Inf,
    ngrid=100,
)
    normal = SVector{3}(normal)
    position = SVector{3}(position)
    if normal ≈ SA_F64[0, 0, 1]
        transform = CoordinateTransformations.Translation(position)
    else
        rot_axis = cross(SA_F64[0, 0, 1], normal)
        rot_angle = acos(normal[3] / norm(normal))
        transform = CoordinateTransformations.AffineMap(
            Rotations.AngleAxis(rot_angle, rot_axis...),
            position
        )
    end
    PlotPlane(CoordinateTransformations.inv(transform), sizex, sizey, ngrid)
end

"""
    transform(plot_plane, position)
    transform(plot_plane, snapshot)

Transforms the `position` (length-3 vector) or [`Snapshot`](@ref) to a space, where the [`PlotPlane`](@ref) lies in the x-y-plane centered on origin.
"""
function transform(pp::PlotPlane, pos::PosVector)
    base = pp.transformation(pos)
    sizex, sizey = pp.sizex, pp.sizey
    correct = [isfinite(sizex) ? sizex : 0., isfinite(sizey) ? sizey : 0., 0.] / 2.
    mod.(base .+ correct, (sizex, sizey, Inf)) .- correct
end

function transform(pp::PlotPlane, snap::Snapshot)
    Snapshot(
        [Spin(transform(pp, s.position), s.orientation) for s in snap.spins],
        snap.time
    )
end


"""
    project_on_grid(plot_plane, snap)

Spins from the [`Snapshot`](@ref) are projected onto the [`PlotPlane`](@ref) in two ways:
- along the normal spins are projected onto the plane from infinitely far (TODO: give finite extent)
- in the other directions any spins are projected onto the plane using mod(position[1], `sizex`) and mod(position[2], `sizey`).
    This assumes that the geometry and field repeats itself ad infinitum beyond the `PlotPlane` (TODO: allow this assumption to be turned off).
In effect, this means that all spins are projected onto the `PlotPlane`.
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


@Makie.recipe(SnapshotPlanarPlot, snap, plane) do scene
    Makie.Attributes(
        markersize=0.1,
        marker=:circle,
        markerspace=:pixel,
        vector=0.1,
        arrowsize=0.1,
    )
end

"""
    plot(snapshot, plot_plane)
    plot!(snapshot, plot_plane)

Plots the spins in the [`Snapshot`](@ref) projected onto given [`PlotPlane`](@ref).
Each spin is represented by an arrow showing the transverse component of the spin.
The background color shows the average spin in the neighbourhood color coded as described in [`color`](@ref).
"""
function Makie.plot!(sp::SnapshotPlanarPlot)
    snap = sp[1]
    planar = sp[2]

    projection = @lift project_on_grid($planar, $snap)
    xs = @lift $projection[1]
    ys = @lift $projection[2]
    c = @lift color.($projection[3])
    Makie.image!(sp, xs, ys, c)

    pos = @lift [Makie.Point2f(s.position[1:2]) for s in transform($planar, $snap).spins]
    vl = sp[:vector]
    directions = @lift [Makie.Point2f(vector(s.orientation)[1:2] .* $vl) for s in $snap.spins]
    kwargs = Dict([sym => sp[sym] for sym in [:arrowsize]] )
    Makie.arrows!(sp, pos, directions; color=:black, kwargs...)
    sp
end

Makie.plottype(::Snapshot, ::PlotPlane) = SnapshotPlanarPlot