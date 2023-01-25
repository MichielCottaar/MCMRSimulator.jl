include("components.jl")
include("gradients.jl")

"""
    Sequence(;TR, gradients=nothing, pulses=nothing, B0=3., interplate_gradients=:step)

An MR sequence represented by a series of pulses repeated with a given repetition time (`TR`).

Possible sequence components are:
- [`InstantRFPulse`](@ref): instantaneous radio-frequency pulse flipping the spin orientations.
- [`InstantGradient`](@ref): instantaneous gradients encoding spatial patterns in the spin phase distribution.
- [`Readout`](@ref): Store the spins at this timepoint.

The index of the next pulse is given by [`next_pulse`](@ref)(sequence, current_time).
The time of this pulse can then be extracted using [`time`](@ref)(sequence, index).
Note that these indices go on till infinite reflecting the repetitive nature of RF pulses over the `TR` time.
"""
struct Sequence{NI, NR, G<:MRGradients}
    scanner :: Scanner
    gradient :: G
    instants :: SVector{NI, InstantComponent}
    TR :: Float
    readout_times :: SVector{NR, Float}
    function Sequence(scanner::Scanner, gradients::MRGradients, pulses::AbstractVector, TR :: Real)
        all([p isa Union{InstantComponent, Readout} for p in pulses])
        instants = sort([p for p in pulses if p isa InstantComponent], by=get_time)
        @assert length(instants) == 0 || (get_time(instants[end]) <= TR && get_time(instants[1]) >= 0)
        readout_times = sort([get_time(p) for p in pulses if p isa Readout])
        @assert length(readout_times) == 0 || (readout_times[end] <= TR && readout_times[1] >= 0)
        result = new{length(instants), length(readout_times), typeof(gradients)}(scanner, gradients, instants, Float(TR), readout_times)
        @assert length(result.instants) == 0 || result.instants[end].time <= TR
        result
    end
end

function Sequence(; scanner=nothing, gradients=[], pulses::AbstractVector=[], TR::Real, B0::Real=3., interpolate_gradients=:linear)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    if !isa(gradients, MRGradients)
        gradients = create_gradients(gradients, TR, interpolate=interpolate_gradients)
    end
    Sequence(scanner, gradients, pulses, TR)
end

Base.getindex(s :: Sequence, index :: Integer) = s.instants[mod(index - 1, length(s.instants)) + 1]
Base.getindex(s :: Sequence{0}, index :: Integer) = error("Can not index an empty sequence")


get_time(sequence :: Sequence{0}, index :: Integer) = Inf
function get_time(sequence :: Sequence, index :: Integer)
    nTR = div(index - 1, length(sequence.instants))
    sequence[index].time + nTR * sequence.TR
end

"""
    next_pulse(sequence, current_time)

Find the index of the next pulse that will be applied.
For an empty sequence 0 will be returned.
"""
next_pulse(sequence :: Sequence{0}, time :: AbstractFloat) = 0
function next_pulse(sequence :: Sequence, time :: AbstractFloat)
    nTR = Int(div(time, sequence.TR, RoundDown))
    within_TR = mod(time, sequence.TR)
    for (index, pulse) in enumerate(sequence.instants)
        if pulse.time >= within_TR
            return index + nTR * length(sequence.instants)
        end
    end
    return 1 + (nTR + 1) * length(sequence.instants)
end

get_gradient(position::AbstractVector, seq::Sequence, time::Number) = get_gradient(position, seq.gradient, mod(time, seq.TR))
get_gradient(position::AbstractVector, seq::Sequence, time1::Number, time2::Number) = get_gradient(position, seq.gradient, mod(time1, seq.TR), mod(time2, seq.TR))
get_gradient(seq::Sequence, time::Number) = get_gradient(seq.gradient, mod(time, seq.TR))
get_gradient(seq::Sequence, time1::Number, time2::Number) = get_gradient(seq.gradient, mod(time1, seq.TR), mod(time2, seq.TR))

include("timestep.jl")
include("diffusion.jl")