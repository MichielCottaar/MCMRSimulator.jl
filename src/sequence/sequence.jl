include("components.jl")

struct Sequence{N, P<:SequenceComponent, T<:AbstractFloat}
    pulses :: SVector{N, P}
    TR :: T
    B0 :: T
    function Sequence(pulses::AbstractVector{<:SequenceComponent}, TR :: Real, B0 :: Real = 3.)
        result = new{length(pulses), eltype(pulses), typeof(TR)}(sort(pulses, by=x->x.time), TR, B0)
        if length(result.pulses) > 0
            @assert result.pulses[end].time <= TR
        end
        result
    end
end

Sequence(TR :: Real, B0 :: Real = 3.) = Sequence(SequenceComponent[], TR, B0)
Base.getindex(s :: Sequence, index :: Integer) = s.pulses[mod(index - 1, length(s.pulses)) + 1]
Base.getindex(s :: Sequence{0}, index :: Integer) = error("Can not index an empty sequence")


time(sequence :: Sequence{0}, index :: Integer) = Inf
function time(sequence :: Sequence, index :: Integer)
    nTR = div(index - 1, length(sequence.pulses))
    sequence[index].time + nTR * sequence.TR
end

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