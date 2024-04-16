module SequenceParts
import StaticArrays: SVector
import LinearAlgebra: norm
import MRIBuilder: BaseSequence, BaseBuildingBlock, waveform_sequence, events, get_gradient, gradient_strength, duration, edge_times, get_pulse, gradient_strength3, slew_rate3, iter_instant_gradients, iter_instant_pulses, TR
import MRIBuilder.Components: NoGradient, ConstantGradient, ChangingGradient, GenericPulse, InstantGradient3D, InstantPulse
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

abstract type NoPulsePart <: SequencePart end

struct EmptyPart <: NoPulsePart
end


struct ConstantPart <: NoPulsePart
    strength :: SVector{3, Float64}
end


struct LinearPart <: NoPulsePart
    start :: SVector{3, Float64}
    final :: SVector{3, Float64}
end


struct PulsePart{T<:NoPulsePart} <: SequencePart
    pulse :: GenericPulse
    gradient :: T
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
        nsteps = isinf(use_timestep) ? 1 : Int(div(t2 - t1, use_timestep, RoundUp))
        append!(splits, range(t1, t2, length=nsteps+1))
    end
    return unique(splits)
end


"""
    MultSequencePart(duration, parts)

A set of N [`SequencePart`](@ref) objects representing overlapping parts of `N` sequences.
"""
struct MultSequencePart{N, T<:SequencePart}
    duration :: Float64
    parts :: SVector{N, T}
    function MultSequencePart(duration::Number, parts::AbstractVector{<:SequencePart})
        N = length(parts)
        if iszero(N)
            return new{0, EmptyPart}(Float64(duration), EmptyPart[])
        end
        base_type = typeof(parts[1])
        if all(p -> p isa base_type, parts)
            T = base_type
        else
            T = SequencePart
        end
        return new{N, T}(
            Float64(duration),
            SVector{N, T}(parts)
        )
    end
end

"""
    SplitSequence(simulation, t1, t2)
    SplitSequence(sequences, t1, t2, timestep)

Splits each sequence between times `t1` and `t2`.
"""
struct SplitSequence{N}
    parts :: Vector{MultSequencePart{N}}
    instants :: Vector{SVector{N, <:Union{Nothing, InstantGradient3D, InstantPulse}}}
end

function SplitSequence(sequences::AbstractVector{<:BaseSequence}, tstart::Number, tfinal::Number, timestep::TimeStep)
    N = length(sequences)
    if iszero(N)
        ts = timestep.max_timestep
        Nsteps = isinf(ts) ? 1 : Int(div(tfinal - tstart, ts, RoundUp))
        duration = (tfinal - tstart) / Nsteps
        msp = MultSequencePart(duration, SequencePart[])
        no_instants = SVector{0, Nothing}()
        return SplitSequence{0}([msp for _ in 1:Nsteps], [no_instants for _ in 1:Nsteps])
    end
    times = split_times(sequences, tstart, tfinal, timestep)
    flipped_parts = [split_into_parts(seq, times) for seq in sequences]
    return SplitSequence{N}(
        [MultSequencePart(times[j+1] - times[j], [flipped_parts[i][j] for i in 1:N]) for j in eachindex(flipped_parts[1])],
        get_instants_array(sequences, times),
    )
end

"""
    split_into_parts(sequence, times)

Splits a sequence into a series of [`SequencePart`](@ref) objects.
"""
function split_into_parts(sequence::BaseSequence{N}, times::AbstractVector{<:Number}) where {N}
    res = SequencePart[]
    for (t1, t2) in zip(times[1:end-1], times[2:end])
        tmean = (t1 + t2) / 2.
        gp = get_pulse(sequence, tmean)
        (gradient, _) = get_gradient(sequence, tmean)
        if gradient isa NoGradient
            grad_part = EmptyPart()
        elseif gradient isa ConstantGradient
            grad_part = ConstantPart(gradient_strength3(gradient))
        elseif gradient isa ChangingGradient
            grad_part = LinearPart(gradient_strength3(sequence, t1), gradient_strength3(sequence, t2))
        else
            error("Gradient waveform $gradient is not implemented in the MCMR simulator yet.")
        end
        if isnothing(gp)
            push!(res, grad_part)
        else
            (pulse, pulse_tmean) = gp
            pulse_t1 = pulse_tmean - tmean + t1
            pulse_t2 = pulse_tmean - tmean + t2
            as_generic = GenericPulse(pulse, pulse_t1, pulse_t2)

            push!(res, PulsePart{typeof(grad_part)}(as_generic, grad_part))
        end
    end
    return res
end

function get_instants_array(sequences::AbstractVector{<:BaseSequence}, times::AbstractVector{<:Number})
    orig_instants = [get_instants_array(seq, times) for seq in sequences]
    flipped_instants = [[orig_instants[i][j] for i in eachindex(sequences)] for j in eachindex(times)]
    return [Union{typeof.(instants)...}[instants...] for instants in flipped_instants]
end

function get_instants_array(sequence::BaseSequence, times::AbstractVector{<:Number})
    res = Any[nothing for _ in times]

    TR1 = div(times[1], TR(sequence), RoundDown)
    TR2 = div(times[end], TR(sequence), RoundUp)
    for (t_instant, instant) in [iter_instant_gradients(sequence)..., iter_instant_pulses(sequence)...]
        for current_TR in TR1:TR2
            t_total = t_instant + current_TR * TR(sequence)
            if (t_total < times[1]) || (t_total > times[end])
                continue
            end
            (_, index) = findmin(t -> abs(t - t_total), times)
            res[index] = instant
        end
    end
    return res
end


end