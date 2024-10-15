module SequenceParts
import StaticArrays: SVector, SizedVector
import LinearAlgebra: norm
import MRIBuilder: BaseSequence, BaseBuildingBlock, waveform_sequence, events, get_gradient, edge_times, get_pulse, iter_instant_gradients, iter_instant_pulses, make_generic, variables
import MRIBuilder.Components: NoGradient, ConstantGradient, ChangingGradient, InstantGradient, InstantPulse, split_timestep
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


"""
    ConstantPulse(amplitude, phase, frequency)

Stores the RF pulse state during a single timestep.

It contains:
- `amplitude` of RF pulse in kHz (assumed to be constant over the timestep).
- `phase` of the RF pulse at the beginning of the timestep in degrees.
- off-resonance `frequency` of the RF pulse in kHz.
"""
struct ConstantPulse
    amplitude :: Float64
    phase :: Float64
    frequency :: Float64
end

struct PulsePart{T<:NoPulsePart} <: SequencePart
    pulse :: Vector{ConstantPulse}
    gradient :: T
end


"""
    MultSequencePart(duration, parts)

A set of N [`SequencePart`](@ref) objects representing overlapping parts of `N` sequences.
"""
struct MultSequencePart{N, T<:SequencePart}
    duration :: Float64
    parts :: SizedVector{N, T}
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
            SizedVector{N, T}(parts)
        )
    end
end

"""
    InstantSequencePart(instants)

A set of `N` instant pulses/gradients that should be applied to the spins.

Some of the instants might be `nothing`.
"""
struct InstantSequencePart{N, T}
    instants :: SizedVector{N, T}
    function InstantSequencePart(instants::AbstractVector)
        T = Union{typeof.(instants)...}
        N = length(instants)
        return new{N, T}(SizedVector{N, T}(instants))
    end
end


function iter_building_blocks_raw(seq::BaseSequence; repeat=false) 
    iter = Iterators.flatten(iter_building_blocks_raw.(seq))
    if repeat
        return Iterators.cycle(iter)
    else
        return iter
    end
end
iter_building_blocks_raw(bb::BaseBuildingBlock) = [bb]
iter_building_blocks(seq::BaseSequence; kwargs...) = Iterators.accumulate(iter_building_blocks_raw(seq; kwargs...); init=(0., nothing)) do (prev_time, prev_bb), bb
    if isnothing(prev_bb)
        return (prev_time, bb)
    else
        return (prev_time + variables.duration(prev_bb), bb)
    end
end
function iter_building_blocks(seq::BaseSequence, tstart, tend)
    after_tstart = Iterators.dropwhile(iter_building_blocks(seq; repeat=true)) do (time, bb)
        time + variables.duration(bb) < tstart
    end
    return Iterators.takewhile(after_tstart) do (time, _)
        time < tend
    end
end

function iter_part_times(seq::BaseSequence, tstart)
    Iterators.flatten(
        Iterators.map(iter_building_blocks(seq, tstart, Inf)) do (time, bb)
            et = edge_times(bb)
            et_with_instant = Tuple{Float64, Any}[(time, nothing) for time in et]
            for (t_instant, instant) in [iter_instant_gradients(bb)..., iter_instant_pulses(bb)...]
                index = findlast(et_with_instant) do (t, _)
                    t == t_instant
                end
                insert!(et_with_instant, index + 1, (t_instant, instant))
            end
            return ((time, t, t2, bb, instant) for ((t, _), (t2, instant)) in zip(et_with_instant[1:end-1], et_with_instant[2:end]))
        end
    )
end

struct _IterEdges{T}
    individual_iterators::Vector{T}
    tstart::Float64
    tfinal::Float64
end

Base.iterate(ie::_IterEdges) = Base.iterate(ie, (ie.tstart, iterate.(ie.individual_iterators)))

function Base.iterate(ie::_IterEdges, my_state)
    current_time, states = my_state

    new_states = map(ie.individual_iterators, states) do iter, state
        (time, _, t2, _, instant), _ = state
        while isnothing(instant) && (time + t2 <= current_time || isapprox(time + t2, current_time, atol=1e-8))
            state = Base.iterate(iter, state[2])
            (time, _, t2, _, instant), _ = state
        end
        return state
    end
    if any(!isnothing(state[1][5]) for state in new_states)
        if current_time ≈ ie.tfinal
            return nothing
        end
        result = (map(new_states) do state
            (_, _, _, _, instant), _ = state
            return instant
        end)
        new_states = map(ie.individual_iterators, new_states) do iter, state
            if isnothing(state[1][5])
                return state
            else
                return Base.iterate(iter, state[2])
            end
        end
        return (result, (current_time, new_states))
    end

    next_time = min(ie.tfinal, 
        minimum(new_states; init=Inf) do ((time, _, t2, _), _)
            time + t2
        end,
    )
    if (next_time ≈ current_time) && (current_time ≈ ie.tfinal)
        return nothing
    else
        return ((current_time, next_time, map(new_states) do state
            (time, _, _, bb), _ = state
            return (time, bb)
        end), (next_time, new_states))
    end
end

function iter_part_times(sequences::AbstractVector{<:BaseSequence}, tstart, tfinal)
    iter_part_times(collect(sequences), tstart, tfinal)
end

function iter_part_times(sequences::Vector{<:BaseSequence}, tstart, tfinal)
    iters = iter_part_times.(sequences, tstart)
    dropped = [Iterators.dropwhile(i) do (time, _, t2, _, _)
        time + t2 < tstart
    end for i in iters]
    return _IterEdges(
        dropped,
        Float64(tstart),
        Float64(tfinal),
    )
end

function iter_part_times(sequences::AbstractVector{<:BaseSequence}, tstart, tfinal, timestep::TimeStep)
    Iterators.flatten(
        Iterators.map(iter_part_times(sequences, tstart, tfinal)) do var
            if !(var isa Tuple)
                return (var, )
            end
            
            (t1, t2, bbs) = var
            tmean = (t1 + t2) / 2

            max_grad = iszero(length(sequences)) ? 0. : maximum(bbs) do (t_bb, bb)
                norm(variables.gradient_strength(get_gradient(bb, tmean - t_bb)[1]))
            end
            use_timestep = timestep(max_grad)
            nsteps = isinf(use_timestep) ? 1 : Int(div(t2 - t1, use_timestep, RoundUp))
            splits = range(t1, t2, length=nsteps+1)
            return (
                (t1p, t2p, bbs)
                for (t1p, t2p) in zip(splits[1:end-1], splits[2:end])
            )
        end
    )
end


function iter_parts(sequences::AbstractVector{<:BaseSequence}, tstart, tfinal, timestep::TimeStep)
    Iterators.map(iter_part_times(sequences, tstart, tfinal, timestep)) do var
        if !(var isa Tuple)
            return InstantSequencePart(var)
        else
            (t1, t2, bbs) = var
            tmean = (t1 + t2) / 2.
            if iszero(length(sequences))
                return MultSequencePart(t2 - t1, EmptyPart[])
            end
            return MultSequencePart(t2 - t1, map(bbs) do (t_bb, bb)
                time_in_bb = tmean - t_bb
                gp = get_pulse(bb, time_in_bb)
                (gradient, _) = get_gradient(bb, time_in_bb)
                if gradient isa NoGradient
                    grad_part = EmptyPart()
                elseif gradient isa ConstantGradient
                    grad_part = ConstantPart(variables.gradient_strength(gradient))
                elseif gradient isa ChangingGradient
                    grad_part = LinearPart(variables.gradient_strength(bb, max(t1 - t_bb, 0.)), variables.gradient_strength(bb, min(t2 - t_bb, variables.duration(bb))))
                else
                    error("Gradient waveform $gradient is not implemented in the MCMR simulator yet.")
                end
                if isnothing(gp)
                    return grad_part
                else
                    (pulse, pulse_tmean) = gp
                    pulse_t1 = pulse_tmean - tmean + t1
                    pulse_t2 = pulse_tmean - tmean + t2
        
                    max_ts = split_timestep(pulse, 1e-3)
        
                    nparts = isinf(max_ts) ? 1 : Int(div(pulse_t2 - pulse_t1, max_ts, RoundUp))
                    ts = (pulse_t2 - pulse_t1) / nparts
        
                    time_parts = ((1:nparts) .- 0.5) .* ts .+ pulse_t1
                    parts = [
                        ConstantPulse(variables.amplitude(pulse, t), variables.phase(pulse, t) - variables.frequency(pulse, t) * 180 * ts, variables.frequency(pulse, t))
                        for t in time_parts
                    ]
                    return PulsePart{typeof(grad_part)}(parts, grad_part)
                end
            end)
        end
    end
end


end