include("components.jl")

struct Sequence{N, P<:SequenceComponent, T<:AbstractFloat}
    pulses :: SVector{N, P}
    TR :: T
    B0 :: T
    function Sequence(pulses::AbstractVector{<:SequenceComponent}, TR :: Real, B0 :: Real = 3.)
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