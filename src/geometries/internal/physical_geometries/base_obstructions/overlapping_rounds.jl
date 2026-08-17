struct OverlappingRound{N} <: BaseObstruction{N}
    radius::Float64
    overlaps_with::Vector{Tuple{SVector{N, Float64}, Float64}}
end

has_inside(::Type{<:OverlappingRound}) = true
has_single_inside(::Type{<:OverlappingRound}) = true

function inside_indices(
    round::OverlappingRound{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    !Base.isempty(intersection) && return intersection.inside ? [ObstructionIndex()] : ObstructionIndex[]
    isinside(round, position) ? [ObstructionIndex()] : ObstructionIndex[]
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
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    candidate = detect_intersection(Round{N}(round.radius), start, destination, previous_hit)
    Base.isempty(candidate) && return candidate
    position = start + candidate.distance .* (destination - start)
    return Intersection(
        candidate.distance,
        candidate.normal,
        candidate.inside,
        ObstructionIndex(),
        _inside_other_rounds(round, position),
    )
end

function surface_sampling(
    round::OverlappingRound{N}, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{N}, scale_density,
) where {N}
    positions, normals = surface_sampling(Round{N}(round.radius), density, bounding_box, scale_density)
    keep = [!_inside_other_rounds(round, position) for position in positions]
    positions[keep], normals[keep]
end

size_scale(round::OverlappingRound) = round.radius
