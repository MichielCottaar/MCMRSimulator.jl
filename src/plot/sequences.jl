module Sequences
using Makie
import Colors
import LinearAlgebra: norm
import ...Sequences: Sequence, MRGradients, RFPulse, Readout, InstantComponent, InstantGradient, InstantRFPulse
import ...Sequences: flip_angle, qval, control_points, gradient
import ...Methods: get_time
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
        max_angle = maximum([abs(flip_angle(p)) for p in s.instants if isa(p, InstantRFPulse)], init=0.)
        max_rf = maximum(Vector([maximum(p.amplitude.amplitudes) for p in s.pulses]), init=0.)
        max_qval = maximum([qval(p) for p in s.instants if isa(p, InstantGradient)], init=0.)
        for pulse in [s.instants..., s.pulses..., Readout.(s.readout_times)...]
            pulseplot!(sp, pulse; max_rf_pulse=max_rf, max_rf_angle=max_angle, max_qval=max_qval)
        end
        max_G = isnothing(max_G) ? s.scanner.gradient : max_G
        if isinf(max_G)
            if length(s.gradients) > 0
                max_G = maximum([maximum(norm.(grad.shape.amplitudes)) for grad in s.gradients])
            end
        end
        for grad in s.gradients
            gradientplot!(sp, grad, max_G=max_G, single_gradient=sg)
        end
    end
    seq[] = seq[]
    sp
end

Makie.plottype(::Sequence) = Sequence_Plot


@Makie.recipe(PulsePlot, pulse) do scene
    Makie.Theme(
        max_rf_pulse=nothing,
        max_rf_angle=180.,
        max_qval=nothing,
    )
end

function Makie.plot!(pp::PulsePlot)
    p = pp[1]
    r = pp[:max_rf_pulse]
    a = pp[:max_rf_angle]
    q = pp[:max_qval]
    comb = @lift ($p, $r, $a, $q)
    on(comb) do as_tuple
        pulse, max_rf, max_angle, max_qval = as_tuple
        if isa(pulse, InstantRFPulse)
            if ~iszero(flip_angle(pulse))
                height = 0.9 * flip_angle(pulse) / max_angle
                Makie.arrows!(pp, [get_time(pulse)], [0.], [0.], [height])
                Makie.text!(pp, "α=" * string(Int(round(flip_angle(pulse)))), position=(get_time(pulse), height + 0.05), align=(:center, :center))
            else
                Makie.text!(pp, string("α=0"), position=(get_time(pulse), 0.05), align=(:center, :center))
            end
        elseif isa(pulse, InstantGradient)
            if ~iszero(qval(pulse))
                height = (isnothing(max_qval) ? 1. : qval(pulse) / max_qval) * 0.9
                Makie.barplot!(pp, [get_time(pulse)], [height])
                Makie.text!(pp, "q=" * string(round(qval(pulse), sigdigits=2)), position=(get_time(pulse), height + 0.05), align=(:center, :center))
            else
                Makie.text!(pp, string("q=0"), position=(get_time(pulse), 0.05), align=(:center, :center))
            end
        elseif isa(pulse, Readout)
            Makie.arrows!(pp, [get_time(pulse)], [0.5], [0.], [-0.5])
            Makie.text!(pp, "readout", position=(get_time(pulse), 0.55), align=(:center, :center))
        elseif isa(pulse, RFPulse)
            times = control_points(pulse.amplitude)
            norm = isnothing(max_rf) ? maximum(pulse.amplitude.amplitudes) : max_rf
            Makie.lines!(pp, times, pulse.amplitude.amplitudes ./ norm, color="black")
        end
    end
    p[] = p[]
    pp
end

Makie.plottype(::InstantComponent) = PulsePlot
Makie.plottype(::Readout) = PulsePlot
Makie.plottype(::RFPulse) = PulsePlot


@Makie.recipe(GradientPlot, pulse) do scene
    Makie.Theme(
        max_G=nothing,
        single_gradient=false,
    )
end

function Makie.plot!(gp::GradientPlot)
    comb = @lift ($(gp[1]), $(gp[:max_G]), $(gp[:single_gradient]))

    on(comb) do as_tuple
        (gradient_profile, max_G, single_grad) = as_tuple
        cp = sort(unique(control_points(gradient_profile)))
        times = sort(unique([
            prevfloat.(cp[2:end])...
            nextfloat.(cp[1:end-1])...
        ]))
        gradients = [gradient(gradient_profile, t) for t in times]
        if length(gradients) > 0
            grad_sizes = [norm(g) for g in gradients]
            if isnothing(max_G) | !isfinite(max_G)
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
    end
    gp[1][] = gp[1][]
    gp
end

Makie.plottype(::MRGradients) = GradientPlot

end