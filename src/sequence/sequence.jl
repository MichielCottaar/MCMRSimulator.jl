include("instants.jl")
include("shape.jl")
include("gradients.jl")
include("radio_frequency.jl")

"""
    Sequence(;TR, gradients=nothing, pulses=nothing, B0=3., interplate_gradients=:step)

An MR sequence represented by a series of pulses repeated with a given repetition time (`TR`).

Possible sequence components are:
- [`InstantRFPulse`](@ref): instantaneous radio-frequency pulse flipping the spin orientations.
- [`InstantGradient`](@ref): instantaneous gradients encoding spatial patterns in the spin phase distribution.
- [`Readout`](@ref): Store the spins at this timepoint.

The index of the next instantaneous pulse is given by [`next_instant`](@ref)(sequence, current_time).
The time of this pulse can then be extracted using [`time`](@ref)(sequence, index).
Note that these indices go on till infinite reflecting the repetitive nature of RF pulses over the `TR` time.
"""
struct Sequence{NI, NP, NR, G<:MRGradients}
    scanner :: Scanner
    gradient :: G
    instants :: SVector{NI, InstantComponent}
    pulses :: SVector{NP, RFPulse}
    TR :: Float
    readout_times :: SVector{NR, Float}
    function Sequence(scanner::Scanner, gradients::MRGradients, pulses::AbstractVector, TR :: Real)
        all([p isa Union{InstantComponent, Readout, RFPulse} for p in pulses])
        instants = sort([p for p in pulses if p isa InstantComponent], by=get_time)
        @assert length(instants) == 0 || (get_time(instants[end]) <= TR && get_time(instants[1]) >= 0)
        @assert length(instants) == length(unique(get_time.(instants)))
        readout_times = sort([get_time(p) for p in pulses if p isa Readout])
        @assert length(readout_times) == 0 || (readout_times[end] <= TR && readout_times[1] >= 0)
        rf_pulses = sort([p for p in pulses if p isa RFPulse], by=start_time)
        for (p1, p2) in zip(rf_pulses[1:end-1], rf_pulses[2:end])
            @assert end_time(p1) <= start_time(p2)
        end
        for pulse in rf_pulses
            @assert start_time(pulse) >= 0
            @assert end_time(pulse) <= TR
        end
        new{length(instants), length(rf_pulses), length(readout_times), typeof(gradients)}(scanner, gradients, instants, rf_pulses, Float(TR), readout_times)
    end
end

function Sequence(; scanner=nothing, gradients=nothing, pulses::AbstractVector=[], TR::Real, B0::Real=3.)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    if isnothing(gradients)
        gradients = MRGradients()
    elseif !isa(gradients, MRGradients)
        gradients = MRGradients(gradients)
    end
    Sequence(scanner, gradients, pulses, TR)
end

function effective_pulse(sequence::Sequence, t0::Number, t1::Number)
    current = current_pulse(sequence, (t0 + t1) / 2)
    if isnothing(current)
        return InstantRFPulse(0, 0, 0)
    else
        @assert t0 >= start_time(current)
        @assert t1 <= end_time(current)
        return effective_pulse(current, t0, t1)
    end
end

start_time(pulse::InstantComponent) = pulse.time
end_time(pulse::InstantComponent) = pulse.time
start_time(number::Number) = number
end_time(number::Number) = number

function add_TR(pulse::InstantComponent, delta_time) 
    @set pulse.time += delta_time
end

next_index(pulses, time) = searchsortedfirst(pulses, time, by=end_time)
previous_index(pulses, time) = searchsortedlast(pulses, time, by=start_time)
function current_index(pulses, time)
    index = next_index(pulses, time)
    return index == previous_index(pulses, time) ? index : -1
end

for relative in ("previous", "current", "next")
    for pulse_type in ("pulse", "instant")
        func_name = Symbol(relative * "_" * pulse_type)
        plural_pulse_type = Symbol(pulse_type * "s")
        index_func = Symbol(relative * "_index")
        @eval begin
            function $func_name(sequence, time)
                nTR = Int(div(time, sequence.TR, RoundDown))
                within_TR = time - sequence.TR * nTR
                all_pulses = sequence.$plural_pulse_type
                if length(all_pulses) == 0
                    return nothing
                end
                index = $index_func(all_pulses, within_TR)
                if index == 0
                    nTR -= 1
                    index = length(all_pulses)
                elseif index == length(all_pulses) + 1
                    nTR += 1
                    index = 1
                end
                if index == -1 || nTR < 0
                    return nothing
                end
                add_TR(all_pulses[index], nTR * sequence.TR)
            end
        end
    end
end

gradient(position::AbstractVector, seq::Sequence, time::Number) = gradient(position, seq.gradient, mod(time, seq.TR))
gradient(position::AbstractVector, seq::Sequence, time1::Number, time2::Number) = gradient(position, seq.gradient, mod(time1, seq.TR), iszero(mod(time2, seq.TR)) ? seq.TR : mod(time2, seq.TR))
gradient(seq::Sequence, time::Number) = gradient(seq.gradient, mod(time, seq.TR))
gradient(seq::Sequence, time1::Number, time2::Number) = gradient(seq.gradient, mod(time1, seq.TR), iszero(mod(time2, seq.TR)) ? seq.TR : mod(time2, seq.TR))

include("timestep.jl")
include("diffusion.jl")