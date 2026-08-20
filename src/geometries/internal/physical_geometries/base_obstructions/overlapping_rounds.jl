struct OverlappingRound{N} <: BaseObstruction{N}
    radius::Float64
    overlaps_with::Vector{Tuple{SVector{N, Float64}, Float64}}
end

has_inside(::Type{<:OverlappingRound}) = true
has_single_inside(::Type{<:OverlappingRound}) = true
Round(r::OverlappingRound{N}) where {N} = Round{N}(r.radius)

distance_to_surface(round::OverlappingRound{N}, position::SVector{N, Float64}) where {N} =
    distance_to_surface(Round(round), position)

function isinside_single(
    round::OverlappingRound{N},
    position::SVector{N, Float64},
    previous_intersection=nothing,
) where {N}
    isinside_single(Round(round), position, previous_intersection)
end

const OverlappingInfiniteCylinder = OverlappingRound{2}
const OverlappingSphere = OverlappingRound{3}

OverlappingRound{N}(radius::Number) where {N} =
    OverlappingRound{N}(Float64(radius), Tuple{SVector{N, Float64}, Float64}[])

InternalBoundingBox(round::OverlappingRound{N}) where {N} = InternalBoundingBox{N}(round.radius)

function _inside_other_rounds(round::OverlappingRound{N}, position::SVector{N, Float64}) where {N}
    any(
        sum((position - offset) .* (position - offset)) < radius^2
        for (offset, radius) in round.overlaps_with
    )
end

function find_intersection(
    round::OverlappingRound{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
) where {N}
    find_intersection(Round{N}(round.radius), start, destination, previous_hit)
end

function get_intersection_params(
    round::OverlappingRound{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    intersection::Tuple,
) where {N}
    params = get_intersection_params(Round{N}(round.radius), start, destination, intersection)
    position = start + intersection[end] .* (destination - start)
    (
        inside=params.inside,
        normal=params.normal,
        hit_gap=_inside_other_rounds(round, position),
    )
end

function surface_sampling(
    round::OverlappingRound{N}, density::GeometryLeafProperties,
    scale_density,
) where {N}
    [
        position
        for position in surface_sampling(Round{N}(round.radius), density, scale_density)
        if !_inside_other_rounds(round, position)
    ]
end

size_scale(round::OverlappingRound) = round.radius
