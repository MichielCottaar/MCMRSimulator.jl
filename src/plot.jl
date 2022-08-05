import Makie: Makie, @lift

color(orient::Union{Spin, SpinOrientation}; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))

@Makie.recipe(SequencePlot, seq) do scene
    Makie.Theme(
    )
end

function Makie.plot!(sp::SequencePlot)
    seq = sp[1]
    times = @lift [p.time for p in $seq.pulses]
    flip_angle = @lift [rad2deg(p.flip_angle) for p in $seq.pulses]
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

function Makie.plot!(sp::SnapshotPlot)
    snap = sp[1]
    colors = @lift color.($snap.spins)
    pos = @lift [Makie.Point3f(s.position) for s in $snap.spins]
    Makie.meshscatter!(sp, pos, color=colors)
    sp
end

Makie.plottype(::Snapshot) = SnapshotPlot


struct PlotPlane
    transformation :: CoordinateTransformations.Transformation
    repeatx :: Real
    repeaty :: Real
end


function PlotPlane(
    normal :: PosVector=SA_F64[0, 0, 1], 
    position :: PosVector=SA_F64[0, 0, 0];
    repeatx::Real=Inf, repeaty::Real=Inf,
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
    PlotPlane(CoordinateTransformations.inv(transform), repeatx, repeaty)
end

function transform(pp::PlotPlane, pos::PosVector)
    mod.(pp.transformation(pos), (pp.repeatx, pp.repeaty, Inf))
end

function transform(pp::PlotPlane, snap::Snapshot)
    Snapshot(
        [Spin(transform(pp, s.position), s.orientation) for s in snap.spins],
        snap.time
    )
end


@Makie.recipe(SnapshotPlanarPlot, snap, plane) do scene
    Makie.Attributes(
        markersize=0.1,
        marker=:circle,
        markerspace=:pixel,
        vector=0.,
        arrowsize=0.1,
    )
end

function Makie.plot!(sp::SnapshotPlanarPlot)
    snap = sp[1]
    planar = sp[2]
    pos = @lift [Makie.Point2f(s.position[1:2]) for s in transform($planar, $snap).spins]
    vl = sp[:vector]
    if iszero(vl[])
        colors = @lift color.($snap.spins)
        kwargs = Dict([sym => sp[sym] for sym in [:markersize, :marker, :markerspace]] )
        Makie.scatter!(sp, pos, color=colors; kwargs...)
    else
        directions = @lift [Makie.Point2f(vector(s.orientation)[1:2] .* $vl) for s in $snap.spins]
        kwargs = Dict([sym => sp[sym] for sym in [:arrowsize]] )
        Makie.arrows!(sp, pos, directions; kwargs...)
    end
    sp
end

Makie.plottype(::Snapshot, ::PlotPlane) = SnapshotPlanarPlot