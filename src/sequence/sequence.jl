include("components.jl")

struct Sequence
    pulses :: Vector{SequenceComponent}
    TR :: Real
    B0 :: Real
    function Sequence(pulses::Vector{T}, TR :: Real, B0 :: Real = 3.) where T <: SequenceComponent
        result = new(sort(pulses, by=x->x.time), TR, B0)
        if length(result.pulses) > 0
            @assert result.pulses[end].time <= TR
        end
        result
    end
end

Sequence(TR :: Real, B0 :: Real = 3.) = Sequence(SequenceComponent[], TR, B0)
Base.getindex(s :: Sequence, index :: Integer) = s.pulses[index]

function time(sequence :: Sequence, index :: Integer)
    if index > length(sequence.pulses)
        return Inf
    end
    return sequence.pulses[index].time
end


include("diffusion.jl")