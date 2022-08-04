import Makie: Makie, @lift

color(orient::Union{Spin, SpinOrientation}; saturation=1., value=1.) = Colors.HSVA(phase(orient) + 180, saturation, value, transverse(orient))

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