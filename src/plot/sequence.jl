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
