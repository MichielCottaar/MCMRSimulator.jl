struct OverlappingRound{N} <: BaseObstruction{N}
    radius::Float64
    overlaps_with::Vector{Tuple{SVector{N, Float64}, Float64}}
end

const OverlappingInfiniteCylinder = OverlappingRound{2}
const OverlappingSphere = OverlappingRound{3}

OverlappingRound{N}(radius::Number) where {N} =
    OverlappingRound{N}(Float64(radius), Tuple{SVector{N, Float64}, Float64}[])

InternalBoundingBox(round::OverlappingRound{N}) where {N} = InternalBoundingBox{N}(round.radius)

isinside(round::OverlappingRound{N}, position::SVector{N, Float64}) where {N} =
    sum(position .* position) < round.radius^2

function _inside_other_rounds(round::OverlappingRound{N}, position::SVector{N, Float64}) where {N}
    any(
        sum((position - offset) .* (position - offset)) < radius^2
        for (offset, radius) in round.overlaps_with
    )
end

function detect_intersection(
    round::OverlappingRound{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    obstruction_index::ObstructionIndex=ObstructionIndex(),
) where {N}
    candidate = detect_intersection(Round{N}(round.radius), start, destination, obstruction_index)
    Base.isempty(candidate) && return candidate
    position = start + candidate.distance .* (destination - start)
    return Intersection(
        candidate.distance,
        candidate.normal,
        candidate.inside,
        candidate.obstruction_index,
        _inside_other_rounds(round, position),
    )
end
