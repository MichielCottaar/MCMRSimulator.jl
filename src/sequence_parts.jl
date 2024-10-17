module SequenceParts
import StaticArrays: SVector, StaticVector, SizedVector
import LinearAlgebra: norm
import MRIBuilder: BaseSequence, BaseBuildingBlock, waveform_sequence, events, get_gradient, edge_times, get_pulse, iter_instant_gradients, iter_instant_pulses, make_generic, variables
import MRIBuilder.Components: GradientWaveform, RFPulseComponent, NoGradient, ConstantGradient, ChangingGradient, InstantGradient, InstantPulse, split_timestep, EventComponent
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

static_vector_type(N) = (N < 50 ? SVector : SizedVector){N}


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
struct MultSequencePart{N, T<:SequencePart, ST<:StaticVector{N, T}}
    duration :: Float64
    parts :: ST
    function MultSequencePart(duration::Number, parts::AbstractVector)
        N = length(parts)
        if iszero(N)
            return new{0, EmptyPart, SVector{0, EmptyPart}}(Float64(duration), SVector{0, EmptyPart}())
        end

        base_type = typeof(parts[1])
        if all(p -> p isa base_type, parts)
            T = base_type
        else
            T = SequencePart
        end
        return new{N, T, static_vector_type(N){T}}(
            Float64(duration),
            static_vector_type(N){T}(parts)
        )
    end
end

"""
    InstantSequencePart(instants)

A set of `N` instant pulses/gradients that should be applied to the spins.

Some of the instants might be `nothing`.
"""
struct InstantSequencePart{N, T, ST<:StaticVector{N, T}}
    instants :: ST
    function InstantSequencePart(instants::AbstractVector)
        T = Union{typeof.(instants)...}
        N = length(instants)
        return new{N, T, static_vector_type(N){T}}(static_vector_type(N){T}(instants))
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
            gradients = Tuple{Float64, Float64, GradientWaveform}[]
            pulses = Tuple{Float64, Float64, RFPulseComponent}[]
            instants = Tuple{Float64, Union{InstantGradient, InstantPulse}}[]

            current_grad = NoGradient(0.)
            current_time = time
            next_time = time
            for key in keys(bb)
                component = bb[key]
                if component isa GradientWaveform
                    duration = component.duration
                    current_time = next_time
                    next_time += duration
                    if duration > 0.
                        push!(gradients, (current_time, next_time, component))
                    end
                elseif component isa Tuple{<:Number, <:EventComponent}
                    delay, event = component
                    if event isa Union{InstantGradient, InstantPulse}
                        push!(instants, (current_time + delay, event))
                    elseif event isa RFPulseComponent
                        duration = variables.duration(event)
                        push!(pulses, (current_time + delay, current_time + delay + duration, event))
                    end
                end
            end
            et = sort!(unique!([
                time,
                [t for (_, t, _) in gradients]...,
                [t for (t, _, _) in pulses]...,
                [t for (_, t, _) in pulses]...,
                [t for (t, _) in instants]...,
            ]))

            current_grad(time) = findfirst(gradients) do (t1, t2, _)
                t2 > time
            end
            current_pulse(time) = findfirst(pulses) do (t1, t2, _)
                t1 <= time < t2
            end

            result = Tuple{Float64, Float64, Tuple{Float64, GradientWaveform}, Union{Nothing, Tuple{Float64, RFPulseComponent}}, Union{Nothing, InstantGradient, InstantPulse}}[]
            for (t1, t2) in zip(et[1:end-1], et[2:end])
                index_grad = current_grad(t1)
                t_grad, _, grad = gradients[index_grad]
                index_pulse = current_pulse(t1)
                pulse = isnothing(index_pulse) ? nothing : (t1 - pulses[index_pulse][1], pulses[index_pulse][end])
                for (t, instant) in instants
                    if t == t1
                        push!(result, (t1, t1, (t1 - t_grad, grad), pulse, instant))
                    end
                end
                push!(result, (t1, t2, (t1 - t_grad, grad), pulse, nothing))
            end
            for (t, instant) in instants
                if t == et[end]
                    t_grad, _, grad = iszero(length(gradients)) ? (t, t, NoGradient(0.)) : gradients[end]
                    push!(result, (t, t, (t - t_grad, grad), nothing, instant))
                end
            end
            return result
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
        (_, t2, _, _, instant), _ = state
        while isnothing(instant) && (t2 <= current_time || isapprox(t2, current_time, atol=1e-8))
            state = Base.iterate(iter, state[2])
            (_, t2, _, _, instant), _ = state
        end
        return state
    end
    if any(!isnothing(state[1][5]) for state in new_states)
        if current_time == ie.tfinal
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
        minimum(new_states; init=Inf) do ((_, t2, _, _, _), _)
            t2
        end,
    )
    if (next_time == current_time) && (current_time == ie.tfinal)
        return nothing
    else
        return ((current_time, next_time, map(new_states) do state
            (t1, _, (t_grad, grad), pulse, _), _ = state
            correct_pulse = isnothing(pulse) ? pulse : (current_time - t1 + pulse[1], pulse[2])
            return ((current_time - t1 + t_grad, grad), correct_pulse)
        end), (next_time, new_states))
    end
end

function iter_part_times(sequences::AbstractVector{<:BaseSequence}, tstart, tfinal)
    iter_part_times(collect(sequences), tstart, tfinal)
end

function iter_part_times(sequences::Vector{<:BaseSequence}, tstart, tfinal)
    iters = iter_part_times.(sequences, tstart)
    dropped = [Iterators.dropwhile(i) do (_, t2, _, _, _)
        t2 < tstart
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
                # instants
                return (var, )
            end
            
            (t1, t2, bbs) = var

            max_grad = iszero(length(sequences)) ? 0. : maximum(bbs) do ((t_grad, grad), _)
                max(
                    norm(variables.gradient_strength(grad, t_grad)),
                    norm(variables.gradient_strength(grad, t_grad + t2 - t1))
                )
            end
            use_timestep = timestep(max_grad)
            nsteps = isinf(use_timestep) ? 1 : Int(div(t2 - t1, use_timestep, RoundUp))
            splits = range(t1, t2, length=nsteps+1)
            return (
                (
                    t1p, t2p, map(bbs) do ((t_grad, grad), pulse)
                        correct_pulse = isnothing(pulse) ? nothing : (pulse[1] + t1p - t1, pulse[2])
                        return ((t_grad + t1p - t1, grad), correct_pulse)
                    end
                )
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
            if iszero(length(sequences))
                return MultSequencePart(t2 - t1, SVector{0, EmptyPart}())
            end
            return MultSequencePart(t2 - t1, map(bbs) do ((t_grad, gradient), gp)
                if gradient isa NoGradient
                    grad_part = EmptyPart()
                elseif gradient isa ConstantGradient
                    grad_part = ConstantPart(variables.gradient_strength(gradient))
                elseif gradient isa ChangingGradient
                    grad_part = LinearPart(variables.gradient_strength(gradient, t_grad), variables.gradient_strength(gradient, t_grad + t2 - t1))
                else
                    error("Gradient waveform $gradient is not implemented in the MCMR simulator yet.")
                end
                if isnothing(gp)
                    return grad_part
                else
                    (pulse_t1, pulse) = gp
                    pulse_t2 = (t2 - t1) + pulse_t1
        
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