include("components.jl")
include("gradients.jl")

"""
    Sequence(;TR, gradients=nothing, pulses=nothing, B0=3., interplate_gradients=:step)

An MR sequence represented by a series of pulses repeated with a given repetition time (`TR`).

Possible sequence components are:
- [`RFPulse`](@ref): instantaneous radio-frequency pulse flipping the spin orientations.
- [`InstantGradient`](@ref): instantaneous gradients encoding spatial patterns in the spin phase distribution.
- [`Readout`](@ref): Store the spins at this timepoint.

The index of the next pulse is given by [`next_pulse`](@ref)(sequence, current_time).
The time of this pulse can then be extracted using [`time`](@ref)(sequence, index).
Note that these indices go on till infinite reflecting the repetitive nature of RF pulses over the `TR` time.
"""
struct Sequence{N, P<:SequenceComponent, G<:MRGradients}
    pulses :: SVector{N, P}
    gradient :: G
    TR :: Float
    B0 :: Float
    function Sequence(gradients::MRGradients, pulses::AbstractVector{<:SequenceComponent}, TR :: Real, B0 :: Real)
        result = new{length(pulses), eltype(pulses), typeof(gradients)}(sort(pulses, by=x->x.time), gradients, Float(TR), Float(B0))
        if length(result.pulses) > 0
            @assert result.pulses[end].time <= TR
        end
        result
    end
end

function Sequence(; gradients=[], pulses::AbstractVector{<:SequenceComponent}=SequenceComponent[], TR::Real, B0::Real=3., interpolate_gradients=:linear)
    if !isa(gradients, MRGradients)
        gradients = create_gradients(gradients, TR, interpolate=interpolate_gradients)
    end
    Sequence(gradients, pulses, TR, B0)
end

Base.getindex(s :: Sequence, index :: Integer) = s.pulses[mod(index - 1, length(s.pulses)) + 1]
Base.getindex(s :: Sequence{0}, index :: Integer) = error("Can not index an empty sequence")


Base.time(sequence :: Sequence{0}, index :: Integer) = Inf
function Base.time(sequence :: Sequence, index :: Integer)
    nTR = div(index - 1, length(sequence.pulses))
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
    for (index, pulse) in enumerate(sequence.pulses)
        if pulse.time >= within_TR
            return index + nTR * length(sequence.pulses)
        end
    end
    return 1 + (nTR + 1) * length(sequence.pulses)
end


include("diffusion.jl")