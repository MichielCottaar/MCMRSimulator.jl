module Sequences
using Makie
import MakieCore
import Colors
import LinearAlgebra: norm
import MCMRSimulator.Sequences: Sequence, MRGradients, RFPulse, Readout, InstantComponent, InstantGradient, InstantRFPulse, read_sequence
import MCMRSimulator.Sequences: flip_angle, qval, control_points, gradient, amplitude, phase
import MCMRSimulator.Methods: get_time
import MCMRSimulator.Plot: Plot_Sequence, print_sequence

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



"""
    plot(sequence)
    plot!(sequence)
    plot_sequence(sequence)
    plot_sequence!(sequence)

Creates a visual representation of a [`Sequence`](@ref) diagram.
"""
function sequence_plot end

function Makie.plot!(scene::Plot_Sequence)
    kwargs = Dict([
        key => scene[key] for key in [
            :visible, :overdraw, :transparency, :fxaa, :inspectable, :depth_shift, :model, :space
        ]
    ])
    text_kwargs = Dict([
        key => scene[key] for key in [
            :textcolor, :font, :fonts, :fontsize
        ]
    ])
    text_kwargs[:color] = map((a, c) -> a === MakieCore.automatic ? c : a, scene[:textcolor], scene[:color])
    line_color = map((a, c) -> a === MakieCore.automatic ? c : a, scene[:linecolor], scene[:color])
    instant_width = map((a, c) -> a * c, scene[:linewidth], scene[:instant_width])

    lift(scene[:sequence]) do sequence
        lines = [
            (label, [0.], [0.], [])
            for label in ("RFx", "RFy", "G1↺", "G2↺", "G3↺", "Gx", "Gy", "Gz")
        ]

        scale = max_finite_pulse(sequence)
        for pulse in sequence.pulses
            times = sort(unique(vcat(pulse.amplitude.times, pulse.phase.times)))
            append!(lines[1][2], times)
            append!(lines[2][2], times)
            append!(lines[1][3], [amplitude(pulse, time) .* cosd.(phase(pulse, time)) ./ scale for time in times])
            append!(lines[2][3], [amplitude(pulse, time) .* sind.(phase(pulse, time)) ./ scale for time in times])
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
            Makie.text!(scene, ln[1] * " ", position=(0., current_y - lower), align=(:right, :center); text_kwargs..., kwargs...)
            Makie.lines!(scene, ln[2], ln[3] .- (lower - current_y); linewidth=scene[:linewidth], color=line_color, kwargs...)
            for (time, height) in ln[4]
                Makie.lines!(scene, [time, time], [current_y - lower, current_y - lower + height]; linewidth=instant_width, color=line_color, kwargs...)
            end

            current_y += (upper - lower) + 0.1
        end

        if length(sequence.readout_times) > 0
            Makie.text!(scene, "ADC "; position=(0., -0.35), align=(:right, :center), text_kwargs..., kwargs...)
            for time in sequence.readout_times
                Makie.lines!(scene, [time, time], [-0.6, -0.1]; color=line_color, linewidth=scene[:linewidth], linestyle=scene[:readout_linestyle], kwargs...)
            end
        end
    end
end

Makie.plottype(::Sequence) = Plot_Sequence

function print_sequence(; sequence_file, output_file, t0, t1, kwargs...)
    sequence = read_sequence(sequence_file)
    f = Figure()
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false)
    plot!(ax, sequence; kwargs...)
    if isinf(t1)
        t1 = sequence.TR
    end
    xlims!(ax, t0 - 0.1 * (t1 - t0), t1)
    ax.xlabel[] = "Time (ms)"
    ax.title[] = sequence_file
    hideydecorations!(ax)
    hidespines!(ax, :l, :r, :t)
    save(output_file, f)
end

end