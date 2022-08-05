import MakieCore: MakieCore, @lift

color(orient::Union{Spin, SpinOrientation}; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))

@MakieCore.recipe(SequencePlot, seq) do scene
    MakieCore.Theme(
    )
end

function MakieCore.plot!(sp::SequencePlot)
    seq = sp[1]
    times = @lift [p.time for p in $seq.pulses]
    flip_angle = @lift [rad2deg(p.flip_angle) for p in $seq.pulses]
    as_zero = @lift $flip_angle * 0
    max_angle = @lift maximum($flip_angle)

    yval = @lift $flip_angle ./ $max_angle
    MakieCore.arrows!(sp, times, as_zero, as_zero, yval)

    text = @lift [string(Int(round(a))) for a in $flip_angle]
    text_positions = @lift [(t, 0.05 + y) for (t, y) in zip($times, $yval)]
    MakieCore.text!(sp, text, position=text_positions, align=(:center, :center))
    #vlines!(sp, TR, color="black")
    sp
end

MakieCore.plottype(::Sequence) = SequencePlot


@MakieCore.recipe(SnapshotPlot, snap) do scene
    MakieCore.Theme(
    )
end

function MakieCore.plot!(sp::SnapshotPlot)
    snap = sp[1]
    colors = @lift color.($snap.spins)
    pos = @lift [MakieCore.Point3f(s.position) for s in $snap.spins]
    MakieCore.meshscatter!(sp, pos, color=colors)
    sp
end

MakieCore.plottype(::Snapshot) = SnapshotPlot


struct PlotPlane
    transformation :: CoordinateTransformations.Transformation
    repeatx :: Real
    repeaty :: Real
    ngrid :: Int
end


function PlotPlane(
    normal :: PosVector=SA_F64[0, 0, 1], 
    position :: PosVector=SA_F64[0, 0, 0];
    repeatx::Real=Inf, repeaty::Real=Inf,
    ngrid=100,
)
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
    PlotPlane(CoordinateTransformations.inv(transform), repeatx, repeaty, ngrid)
end

function transform(pp::PlotPlane, pos::PosVector)
    base = pp.transformation(pos)
    repeatx, repeaty = pp.repeatx, pp.repeaty
    correct = [isfinite(repeatx) ? repeatx : 0., isfinite(repeaty) ? repeaty : 0., 0.] / 2.
    mod.(base .+ correct, (repeatx, repeaty, Inf)) .- correct
end

function transform(pp::PlotPlane, snap::Snapshot)
    Snapshot(
        [Spin(transform(pp, s.position), s.orientation) for s in snap.spins],
        snap.time
    )
end

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


@MakieCore.recipe(SnapshotPlanarPlot, snap, plane) do scene
    MakieCore.Attributes(
        markersize=0.1,
        marker=:circle,
        markerspace=:pixel,
        vector=0.1,
        arrowsize=0.1,
    )
end

function MakieCore.plot!(sp::SnapshotPlanarPlot)
    snap = sp[1]
    planar = sp[2]

    projection = @lift project_on_grid($planar, $snap)
    xs = @lift $projection[1]
    ys = @lift $projection[2]
    c = @lift color.($projection[3])
    MakieCore.image!(sp, xs, ys, c)

    pos = @lift [MakieCore.Point2f(s.position[1:2]) for s in transform($planar, $snap).spins]
    vl = sp[:vector]
    directions = @lift [MakieCore.Point2f(vector(s.orientation)[1:2] .* $vl) for s in $snap.spins]
    kwargs = Dict([sym => sp[sym] for sym in [:arrowsize]] )
    MakieCore.arrows!(sp, pos, directions; color=:black, kwargs...)
    sp
end

MakieCore.plottype(::Snapshot, ::PlotPlane) = SnapshotPlanarPlot