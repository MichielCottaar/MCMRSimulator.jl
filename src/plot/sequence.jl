@Makie.recipe(Sequence_Plot, seq) do scene
    Makie.Theme(
        max_G=nothing,
        single_gradient=false,
    )
end

"""
    plot(sequence)
    plot!(sequence)
    plot_sequence(sequence)
    plot_sequence!(sequence)

Creates a visual representation of a [`Sequence`](@ref) diagram.
"""
function sequence_plot end

function Makie.plot!(sp::Sequence_Plot)
    seq = sp[1]
    on(@lift ($(sp[1]), $(sp[:max_G]), $(sp[:single_gradient]))) do as_tuple
        (s, max_G, sg) = as_tuple
        if any(p->isa(p, RFPulse), s.pulses)
            max_angle = maximum([flip_angle(p) for p in s.pulses if isa(p, RFPulse)])
        else
            max_angle = nothing
        end
        if any(p->isa(p, InstantGradient), s.pulses)
            max_qval = maximum([qval(p) for p in s.pulses if isa(p, InstantGradient)])
        else
            max_qval = nothing
        end
        for pulse in s.pulses
            pulseplot!(sp, pulse; max_rf_pulse=max_angle, max_qval=max_qval)
        end
        gradientplot!(sp, s.gradient, max_G=max_G, single_gradient=sg)
    end
    seq[] = seq[]
    sp
end

Makie.plottype(::Sequence) = Sequence_Plot


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
            Makie.arrows!(pp, [time(pulse)], [0.], [0.], [height])
            Makie.text!(pp, string(Int(round(flip_angle(pulse)))), position=(time(pulse), height + 0.05), align=(:center, :center))
        elseif isa(pulse, InstantGradient)
            height = (isnothing(max_qval) ? 1. : qval(pulse) / max_qval) * 0.9
            Makie.barplot!(pp, [time(pulse)], [height])
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


@Makie.recipe(GradientPlot, pulse) do scene
    Makie.Theme(
        max_G=nothing,
        single_gradient=false,
    )
end

function Makie.plot!(gp::GradientPlot)
    comb = @lift ($(gp[1]), $(gp[:max_G]), $(gp[:single_gradient]))

    on(comb) do as_tuple
        (gradient, max_G, single_grad) = as_tuple
        times = Float[]
        push!(times, nextfloat(Float(gradient.times[1])))
        for t in gradient.times[2:end-1]
            if t > times[end]
                push!(times, prevfloat(Float(t)))
                push!(times, nextfloat(Float(t)))
            end
        end
        push!(times, prevfloat(Float(gradient.times[end])))
        gradients = [get_gradient(gradient, t) for t in times]
        grad_sizes = [norm(g) for g in gradients]
        if isnothing(max_G)
            max_G = maximum(grad_sizes)
        end
        if single_grad
            rgb = [Colors.RGB(abs.(g./s)...) for (g, s) in zip(gradients, grad_sizes)]
            Makie.lines!(gp, times, [s / max_G for s in grad_sizes], color=1:length(times), colormap=rgb)
        else
            for (dim, color) in zip(1:3, ("red", "green", "blue"))
                Makie.lines!(gp, times, [g[dim] / max_G for g in gradients], color=color)
            end
        end
    end
    gp[1][] = gp[1][]
    gp
end

Makie.plottype(::MRGradients) = GradientPlot