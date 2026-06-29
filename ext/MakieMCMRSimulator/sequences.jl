module Sequences
using Makie
import MCMRSimulator.Plot: plot_sequence, plot_sequence!, SequenceDiagram, normalise, SingleSequenceDiagramLine, duration_diagram
import MCMRSimulator.SequenceParts: SequenceWaveform
import MCMRSimulator.Pulseq: PulseqSequence


const SequenceLike = Union{SequenceWaveform, PulseqSequence, SequenceDiagram}

"""
    plot(sequence; kwargs...)
    plot!([scene,] sequence; kwargs...)
    plot_sequence(sequence; kwargs...)
    plot_sequence!([scene,] sequence; kwargs...)

Plots a given sequence in a sequence diagram.

The sequence diagram is a 2D plot with time on the x-axis and the different sequence components on the y-axis. 

This function will only work if a [`Makie`](https://makie.org) backend is imported.
"""
@recipe Plot_Sequence (sequence::SequenceLike, ) begin
    color = @inherit textcolor
    linecolor = Makie.automatic
    linewidth = 1.5
    instant_width = 3.
    textcolor = Makie.automatic
    font = @inherit font
    fonts = @inherit fonts
    fontsize = @inherit fontsize
    Makie.mixin_generic_plot_attributes()...
end

argument_names(::Type{<:Plot_Sequence}, N) = (:sequence, )

Makie.plottype(::SequenceLike) = Plot_Sequence

function Makie.plot!(scene:: Plot_Sequence)

    map!(scene.attributes, [:linewidth, :instant_width], :linewidth_instant) do l, i
        l * i
    end

    map!(scene.attributes, [:textcolor, :color], :textcolor_final) do t, c
        t === Makie.automatic ? c : t
    end

    map!(scene.attributes, [:linecolor, :color], :linecolor_final) do l, c
        l === Makie.automatic ? c : l
    end

    map!(scene.attributes, [:sequence], :sequence_diagram) do sequence
        SequenceDiagram(sequence)
    end

    Makie.register_computation!(scene.attributes, [:sequence_diagram], [:text, :text_x, :text_y, :times, :amplitudes, :event_times, :event_amplitudes, :fake_x, :fake_y]) do inputs, changed, cached
        sequence_diagram = normalise(inputs.sequence_diagram)

        text = String[]
        text_x = Float64[]
        text_y = Float64[]
        times = Float64[]
        amplitudes = Float64[]

        event_times = Float64[]
        event_amplitudes = Float64[]

        current_y = 0.
        for label in (:ADC, :Gz, :Gy, :Gx, :RFy, :RFx)
            line = getproperty(sequence_diagram, label)
            (lower, upper) = range_full(line)
            if lower ≈ upper
                continue
            end
            shift = current_y - lower
            push!(text, string(label))
            push!(text_x, -duration_diagram(sequence_diagram) / 40)
            push!(text_y, shift)
            append!(times, line.times)
            append!(amplitudes, line.amplitudes .+ shift)
            push!(times, NaN)
            push!(amplitudes, NaN)
            for (time, amplitude) in zip(line.event_times, line.event_amplitudes)
                append!(event_times, [time, time, NaN])
                append!(event_amplitudes, [shift, amplitude + shift, NaN])
            end
            current_y += (upper - lower) + 0.1
        end

        fake_x = [-duration_diagram(sequence_diagram) / 10., 0.]
        fake_y = [current_y/2, current_y/2.]

        return (text, text_x, text_y, times, amplitudes, event_times, event_amplitudes, fake_x, fake_y)
    end

    Makie.text!(scene, scene.attributes, scene[:text_x], scene[:text_y]; align=(:right, :center), color=scene[:textcolor_final])
    Makie.lines!(scene, scene.attributes, scene[:times], scene[:amplitudes]; color=scene[:linecolor_final])
    Makie.lines!(scene, scene.attributes, scene[:event_times], scene[:event_amplitudes]; linewidth=scene[:linewidth_instant], color=scene[:linecolor_final])

    Makie.lines!(scene, scene[:fake_x], scene[:fake_y], color=(:black, 0.))
end

range_line(spl::SingleSequenceDiagramLine) = (
    minimum(spl.amplitudes; init=0.),
    maximum(spl.amplitudes; init=0.),
)

range_event(spl::SingleSequenceDiagramLine) = (
    minimum(spl.event_amplitudes; init=0.),
    maximum(spl.event_amplitudes; init=0.),
)

function range_full(spl::SingleSequenceDiagramLine)
    (l1, u1) = range_line(spl)
    (l2, u2) = range_event(spl)
    return (
        min(l1, l2),
        max(u1, u2),
    )
end

end