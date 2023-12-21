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


function Plot.plot_sequence(sequence::Sequence; figure=Dict{Symbol, Any}(), axis=Dict{Symbol, Any}(), kwargs...)
    f = Figure(; figure...)
    ax = Axis(f[1, 1]; axis...)
    plot_sequence!(ax, sequence; kwargs...)
    if length(ax.scene.plots) == 0
        return Makie.FigureAxis(f, ax)
    else
        return Makie.FigureAxisPlot(f, ax, ax.scene.plots[end])
    end
end

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

end