module Sequences
using Makie
import Colors
import LinearAlgebra: norm
import MCMRSimulator.Sequences: Sequence, MRGradients, RFPulse, Readout, InstantComponent, InstantGradient, InstantRFPulse
import MCMRSimulator.Sequences: flip_angle, qval, control_points, gradient
import MCMRSimulator.Methods: get_time
import MCMRSimulator.Plot: Plot, plot_sequence!



max_finite_gradient(sequence::Sequence{L, K, M, 0}) where {L, K, M} = 0.
max_finite_gradient(sequence::Sequence) = maximum(max_finite_gradient.(sequence.gradients))
max_finite_gradient(grad::MRGradients) = maximum(abs.(reduce(vcat, grad.shape.amplitudes)), init=0.)

max_instant_gradient(sequence::Sequence{0}) = 0.
max_instant_gradient(sequence::Sequence) = maximum(max_instant_gradient.(sequence.instants))
max_instant_gradient(grad::InstantGradient) = maximum(abs.(grad.qvec))
max_instant_gradient(pulse::InstantRFPulse) = 0.

max_finite_pulse(sequence::Sequence{L, 0}) where {L} = 0.
max_finite_pulse(sequence::Sequence) = maximum(max_finite_pulse.(sequence.pulses))
max_finite_pulse(pulse::RFPulse) = maximum(abs.(pulse.amplitude.amplitudes), init=0.)

max_instant_pulse(sequence::Sequence{0}) = 0.
max_instant_pulse(sequence::Sequence) = maximum(max_instant_pulse.(sequence.instants))
max_instant_pulse(grad::InstantGradient) = 0.
max_instant_pulse(pulse::InstantRFPulse) = pulse.flip_angle


function plot_sequence!(scene, sequence::Sequence)
    lines = [
        (label, [0.], [0.], [])
        for label in ("RFx", "RFy", "G1↺", "G2↺", "G3↺", "Gx", "Gy", "Gz")
    ]

    scale = max_finite_pulse(sequence)
    for pulse in sequence.pulses
        append!(lines[1][2], pulse.amplitude.times)
        append!(lines[2][2], pulse.amplitude.times)
        append!(lines[1][3], pulse.amplitude.amplitudes .* cosd.(pulse.phase.amplitudes) ./ scale)
        append!(lines[2][3], pulse.amplitude.amplitudes .* sind.(pulse.phase.amplitudes) ./ scale)
    end

    scale = max_instant_pulse(sequence)
    for instant_pulse in sequence.instants
        if !(instant_pulse isa InstantRFPulse)
            continue
        end
        if cosd(instant_pulse.phase) != 0.
            push!(lines[1][4], (instant_pulse.time, instant_pulse.flip_angle * cosd(instant_pulse.phase) / scale))
        end
        if sind(instant_pulse.phase) != 0.
            push!(lines[2][4], (instant_pulse.time, instant_pulse.flip_angle * sind(instant_pulse.phase) / scale))
        end
    end

    instant_scale = max_instant_gradient(sequence)
    finite_scale = max_finite_gradient(sequence)
    for base_index in (2, 5)
        for dim in (1, 2, 3)
            index = base_index + dim

            for grad in sequence.gradients
                if xor(base_index == 2, grad.apply_bvec)
                    continue
                end
                append!(lines[index][2], grad.shape.times)
                append!(lines[index][3], [a[dim] for a in grad.shape.amplitudes] ./ finite_scale)
            end
            for instant_grad in sequence.instants
                if !(instant_grad isa InstantGradient)
                    continue
                end
                if xor(base_index == 2, instant_grad.apply_bvec)
                    continue
                end
                push!(lines[index][4], (instant_grad.time, instant_grad.qvec[dim] ./ instant_scale))
            end
        end
    end

    current_y = 0.
    for ln in lines[end:-1:1]
        push!(ln[2], sequence.TR)
        push!(ln[3], 0.)
        lower = min(minimum(ln[3], init=0.), minimum([a[2] for a in ln[4]], init=0.))
        upper = max(maximum(ln[3], init=0.), maximum([a[2] for a in ln[4]], init=0.))
        if lower == upper
            continue
        end
        Makie.text!(scene, ln[1] * " ", position=(0., current_y - lower), align=(:right, :center))
        Makie.lines!(scene, ln[2], ln[3] .- (lower - current_y), color=:black, linewidth=1.5)
        for (time, height) in ln[4]
            Makie.lines!(scene, [time, time], [current_y - lower, current_y - lower + height], color=:black, linewidth=5)
        end

        current_y += (upper - lower) + 0.1
    end
end

function plot_sequence(sequence::Sequence)
    sc = Scene()
    plot_sequence!(sc, sequence)
    return sc
end


function old()
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
                Makie.arrows!(pp, [get_time(pulse)], [0.], [0.], [height], color=:black)
                Makie.text!(pp, "α=" * string(Int(round(flip_angle(pulse)))), position=(get_time(pulse), height + 0.05), align=(:center, :center))
            else
                Makie.text!(pp, string("α=0"), position=(get_time(pulse), 0.05), align=(:center, :center))
            end
        elseif isa(pulse, InstantGradient)
            if ~iszero(qval(pulse))
                height = (isnothing(max_qval) ? 1. : qval(pulse) / max_qval) * 0.9
                Makie.barplot!(pp, [get_time(pulse)], [height], color=:black)
                Makie.text!(pp, "q=" * string(round(qval(pulse), sigdigits=2)), position=(get_time(pulse), height + 0.05), align=(:center, :center))
            else
                Makie.text!(pp, string("q=0"), position=(get_time(pulse), 0.05), align=(:center, :center))
            end
        elseif isa(pulse, Readout)
            Makie.arrows!(pp, [get_time(pulse)], [0.5], [0.], [-0.5], color=:black)
            Makie.text!(pp, "readout", position=(get_time(pulse), 0.55), align=(:center, :center))
        elseif isa(pulse, RFPulse)
            times = control_points(pulse.amplitude)
            norm = isnothing(max_rf) ? maximum(pulse.amplitude.amplitudes) : max_rf
            Makie.lines!(pp, times, pulse.amplitude.amplitudes ./ norm, color=:black)
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

"""
    plot_gradient!(axis, sequence, dimension)
"""
function plot_gradient!(axis, sequence::Sequence, dimension::Int)
    times = [0.]
    amplitudes = [0.]

    plot_rotation = iszero(dimension)

    for grad in sequence.gradients
        if xor(plot_rotation, ~grad.apply_bvec)
            continue
        end
        shape = grad.shape
        for (time, full_amplitude) in zip(shape.times, shape.amplitudes)
            push!(times, time)
            if plot_rotation
                push!(amplitudes, norm(full_amplitude))
            else
                push!(amplitudes, full_amplitude[dimension])
            end
        end
    end
    push!(times, sequence.TR)
    push!(amplitudes, 0.)
    lines!(axis, times, amplitudes)
    for grad in sequence.instants
        if ~(grad isa InstantGradient)
            continue
        end

    end
end

Makie.plottype(::MRGradients) = GradientPlot

end