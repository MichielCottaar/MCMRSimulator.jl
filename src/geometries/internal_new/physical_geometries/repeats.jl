"""Periodic repetition of a physical geometry."""
module Repeats

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox

struct Repeat{N, P<:PhysicalGeometry{N}} <: PhysicalGeometry{N}
    geometry::P
    repeats::SVector{N, Float64}

    function Repeat{N, P}(geometry::P, repeats::SVector{N, Float64}) where {N, P<:PhysicalGeometry{N}}
        all(repeats .> 0) || throw(ArgumentError("repeat distances must be positive"))
        box = InternalBoundingBox(geometry)
        all(InternalBoundingBoxes.lower(box) .>= -1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far below its repeat bounds"))
        all(InternalBoundingBoxes.upper(box) .<= 1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far above its repeat bounds"))
        new{N, P}(geometry, repeats)
    end
end

Repeat(geometry::P, repeats::AbstractVector{<:Real}) where {N, P<:PhysicalGeometry{N}} =
    Repeat{N, P}(geometry, SVector{N, Float64}(repeats))

has_inside(::Type{<:Repeat{N, P}}) where {N, P} = has_inside(P)

_wrap(repeat::Repeat, position) = mod.(position .+ repeat.repeats / 2, repeat.repeats) .- repeat.repeats / 2

inside_indices(repeat::Repeat{N}, position::SVector{N, Float64}) where {N} =
    inside_indices(repeat.geometry, _wrap(repeat, position))

InternalBoundingBox(::Repeat) = throw(ArgumentError("repeated geometries do not have a finite bounding box"))

function detect_intersection(
    repeat::Repeat{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    displacement = destination - start
    iszero(displacement) && return Intersection{N}()
    closest = Intersection{N}()
    previous = previous_hit
    scaled_start = (start .+ repeat.repeats / 2) ./ repeat.repeats
    scaled_destination = (destination .+ repeat.repeats / 2) ./ repeat.repeats
    for (_, entry_time, _, exit_time, _) in ray_grid_intersections(scaled_start, scaled_destination)
        segment_start = start + entry_time .* displacement
        segment_destination = start + exit_time .* displacement
        local_start = _wrap(repeat, segment_start)
        local_destination = local_start + (segment_destination - segment_start)
        local_intersection = detect_intersection(
            repeat.geometry,
            local_start,
            local_destination,
            previous,
        )
        previous = Intersection{3}()
        Base.isempty(local_intersection) && continue
        distance = entry_time + local_intersection.distance * (exit_time - entry_time)
        if distance < closest.distance
            closest = Intersection(
                distance,
                local_intersection.normal,
                local_intersection.inside,
                local_intersection.obstruction_index,
                local_intersection.hit_gap,
            )
        end
    end
    closest
end

end
