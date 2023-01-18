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
struct Sequence{N, P<:SequenceComponent, G<:MRGradients, M}
    scanner :: Scanner
    gradient :: G
    pulses :: SVector{N, P}
    TR :: Float
    readout_times :: SVector{M, Float}
    function Sequence(scanner::Scanner, gradients::MRGradients, pulses::AbstractVector{<:SequenceComponent}, TR :: Real)
        sorted = sort(pulses, by=x->x.time)
        pulses = [p for p in sorted if !isa(p, Readout)]
        readout_times = [p.time for p in sorted if isa(p, Readout)]
        result = new{length(pulses), eltype(pulses), typeof(gradients), length(readout_times)}(scanner, gradients, pulses, Float(TR), readout_times)
        @assert length(result.pulses) == 0 || result.pulses[end].time <= TR
        result
    end
end

function Sequence(; scanner=nothing, gradients=[], pulses::AbstractVector{<:SequenceComponent}=SequenceComponent[], TR::Real, B0::Real=3., interpolate_gradients=:linear)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    if !isa(gradients, MRGradients)
        gradients = create_gradients(gradients, TR, interpolate=interpolate_gradients)
    end
    Sequence(scanner, gradients, pulses, TR)
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

get_gradient(position, seq::Sequence, time) = get_gradient(position, seq.gradient, mod(time, seq.TR))
get_gradient(position, seq::Sequence, time1, time2) = get_gradient(position, seq.gradient, mod(time1, seq.TR), mod(time2, seq.TR))
get_gradient(seq::Sequence, time) = get_gradient(seq.gradient, mod(time, seq.TR))
get_gradient(seq::Sequence, time1, time2) = get_gradient(seq.gradient, mod(time1, seq.TR), mod(time2, seq.TR))

include("timestep.jl")
include("diffusion.jl")