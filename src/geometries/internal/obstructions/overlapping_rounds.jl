module OverlappingRounds
import StaticArrays: SVector
import Distributions: Poisson
import Random: rand
import ..Rounds: Rounds
import ..ObstructionIntersections: empty_obstruction_intersections, ObstructionIntersection

struct OverlappingRound{N} <: Rounds.FixedObstruction{N}
    radius :: Float64
    overlaps_with :: Vector{Tuple{SVector{N, Float64}, Float64}} # vector of (offset, radius) tuples for other overlapping rounds
end

const OverlappingCylinder = OverlappingRound{2}
const OverlappingSphere = OverlappingRound{3}

OverlappingRound{N}(radius::Number) where {N} = OverlappingRound{N}(radius, Tuple{SVector{N, Float64}, Float64}[])

Rounds.radius(r::OverlappingRound) = r.radius
Rounds.Round{N}(r::OverlappingRound{N}) where {N} = Rounds.Round{N}(r.radius)
Rounds.isinside(r::OverlappingRound{N}, pos::SVector{N}) where {N} = sum(pos .* pos) < Rounds.radius(r)^2
Rounds.has_inside(::Type{<:OverlappingRound}) = true

function inside_other_rounds(overlapping_round::OverlappingRound{N}, position::SVector{N, Float64}) where {N}
    for (offset, radius) in overlapping_round.overlaps_with
        if sum((position - offset) .* (position - offset)) < radius^2
            return true
        end
    end
    return false
end

function Rounds.detect_intersection(r::OverlappingRound{N}, start::SVector{N}, dest::SVector{N}, args...) where {N}
    candidate = Rounds.detect_intersection(Rounds.Round{N}(r), start, dest, args...)
    if candidate === empty_obstruction_intersections[N]
        return candidate
    end
    intersect_location = start + candidate.distance * (dest - start)
    return ObstructionIntersection(
        candidate.distance,
        candidate.normal,
        candidate.inside,
        inside_other_rounds(r, intersect_location)
    )
end


function Rounds.random_surface_positions(r::OverlappingRound{N}, density::Number) where {N}
    positions, normals = Rounds.random_surface_positions(Rounds.Round{N}(r), density)
    # filter out positions that are inside other overlapping rounds
    keep = [!inside_other_rounds(r, pos) for pos in positions]
    return (positions[keep], normals[keep])
end





end