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
    on(seq) do s
        max_angle = maximum([flip_angle(p) for p in s.pulses if isa(p, RFPulse)])
        max_qval = maximum([qval(p) for p in s.pulses if isa(p, InstantGradient)])
        for pulse in s.pulses
            pulseplot!(sp, pulse; max_rf_pulse=max_angle, max_qval=max_qval)
        end
    end
    seq[] = seq[]
    sp
end

Makie.plottype(::Sequence) = SequencePlot


@Makie.recipe(PulsePlot, pulse) do scene
    Makie.Theme(
        max_rf_pulse=180.,
        max_qval=nothing,
    )
end

function Makie.plot!(pp::PulsePlot)
    p = pp[1]
    r = pp[:max_rf_pulse]
    q = pp[:max_qval]
    comb = @lift ($p, $r, $q)
    on(comb) do as_tuple
        pulse, max_rf, max_qval = as_tuple
        if isa(pulse, RFPulse)
            height = 0.9 * flip_angle(pulse) / max_rf
            Makie.arrows!(pp, [time(pulse)], [0.], [0.], [height], markerspace=:pixel, arrowsize=2.)
            Makie.text!(pp, string(Int(round(flip_angle(pulse)))), position=(time(pulse), height + 0.05), align=(:center, :center))
        elseif isa(pulse, InstantGradient)
            height = (isnothing(max_qval) ? 1. : qval(pulse) / max_qval) * 0.9
            Makie.barplot!(pp, time(pulse), height)
            Makie.text!(pp, string(round(qval(pulse), sigdigits=2)), position=(time(pulse), height + 0.05), align=(:center, :center))
        elseif isa(pulse, Readout)
            Makie.arrows!(pp, [time(pulse)], [0.5], [0.], [-0.5])
            Makie.text!(pp, "readout", position=(time(pulse), 0.55), align=(:center, :center))
        end
    end
    p[] = p[]
    pp
end

Makie.plottype(::SequenceComponent) = PulsePlot