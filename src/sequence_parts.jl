module SequenceParts
import StaticArrays: SVector
import LinearAlgebra: norm
import MRIBuilder: BaseSequence, BaseBuildingBlock, waveform_sequence, events, get_gradient, gradient_strength, duration, edge_times, get_pulse, gradient_strength3, slew_rate3
import MRIBuilder.Components: NoGradient, ConstantGradient, ChangingGradient, GenericPulse
import ..TimeSteps: TimeStep


"""
A short part of the sequence that can be handled by the simulator.

There are 4 types:
- [`EmptyPart`](@ref): no pulse or gradient
- [`ContantPart`](@ref): no pulse; constant gradient
- [`LinearPart`](@ref): no pulse; linear gradient
- [`PulsePart`](@ref): RF pulse; constant gradient

To break down a generic sequence in these parts will require some approximations.
"""
abstract type SequencePart end

struct EmptyPart <: SequencePart
    duration :: Float64
end


struct ConstantPart <: SequencePart
    duration :: Float64
    strength :: SVector{3, Float64}
end


struct LinearPart <: SequencePart
    duration :: Float64
    start :: SVector{3, Float64}
    final :: SVector{3, Float64}
end


struct PulsePart <: SequencePart
    pulse :: GenericPulse
    gradient :: SVector{3, Float64}
end

"""
    split_times(sequence(s), timestep)
    split_times(sequence(s), t1, t2, timestep)

Suggests at what times to split one or more sequence into individual timestep for simulation.

The split times will include any time when (for any of the provided sequences):
- the starting/end points of building blocks, gradients, RF pulses, or ADC readouts.
- a gradient or RF pulse is discontinuous in the first derivative
- the time of any instantaneous gradients, RF pulses, or readouts.

Continuous gradient waveforms or RF pulses might be split up further to ensure that the maximum timestep is obeyed.
"""
split_times(sequence::BaseSequence, args...; kwargs...) = split_times([sequence], args...; kwargs...)
split_times(sequences::AbstractVector{<:BaseSequence}, timestep::TimeStep; kwargs...) = split_times(sequences, 0., maximum(duration.(sequences)), timestep; kwargs...)

function split_times(sequences::AbstractVector{<:BaseSequence}, tstart::Number, tfinal::Number, timestep::TimeStep)
    edges = Float64.([tstart, tfinal])
    for sequence in sequences
        raw_edges = edge_times(sequence)
        nTR_start = Int(div(tstart, duration(sequence), RoundDown))
        nTR_final = Int(div(tfinal, duration(sequence), RoundUp))
        for nTR in nTR_start:nTR_final
            for time in raw_edges .+ (nTR * duration(sequence))
                if tstart < time < tfinal
                    push!(edges, time)
                end
            end
        end
    end
    edges = sort(unique(edges))
    splits = Float64[]

    for (t1, t2) in zip(edges[1:end-1], edges[2:end])
        tmean = (t1 + t2) / 2

        max_grad = maximum(map(seq -> norm(gradient_strength(get_gradient(seq, tmean)[1])), sequences))

        use_timestep = timestep(max_grad)
        @show use_timestep
        nsteps = isinf(use_timestep) ? 1 : Int(div(t2 - t1, use_timestep, RoundUp))
        append!(splits, range(t1, t2, length=nsteps+1))
    end
    return unique(splits)
end


"""
    split_into_parts(simulation)
    split_into_parts(sequences, t1, t2, timestep)
    split_into_parts(sequence, times)

Splits a sequence into a series of [`SequencePart`](@ref) objects.
"""
function split_into_parts(sequences::AbstractVector{<:BaseSequence}, tstart::Number, tfinal::Number, timestep::TimeStep)
    times = split_times(sequences, tstart, tfinal, timestep)
    return hcat([split_into_parts(seq, times) for seq in sequences]...)
end

function split_into_parts(sequence::BaseSequence{N}, times::AbstractVector{<:Number}) where {N}
    res = SequencePart[]
    for (t1, t2) in zip(times[1:end-1], times[2:end])
        tmean = (t1 + t2) / 2.
        gp = get_pulse(sequence, tmean)
        (gradient, _) = get_gradient(sequence, tmean)
        if isnothing(gp)
            if gradient isa NoGradient
                push!(res, EmptyPart(t2 - t1))
            elseif gradient isa ConstantGradient
                push!(res, ConstantPart(t2 - t1, gradient_strength3(gradient)))
            elseif gradient isa ChangingGradient
                push!(res, LinearPart(t2 - t1, gradient_strength3(sequence, t1), gradient_strength3(sequence, t2)))
            else
                error("Gradient waveform $gradient is not implemented in the MCMR simulator yet.")
            end
        else
            @assert gradient isa Union{NoGradient, ConstantGradient}
            (pulse, pulse_tmean) = gp
            pulse_t1 = pulse_tmean - tmean + t1
            pulse_t2 = pulse_tmean - tmean + t2
            as_generic = GenericPulse(pulse, pulse_t1, pulse_t2)

            push!(res, PulsePart(as_generic, gradient_strength3(sequence, tmean)))
        end
    end
    return res
end



end