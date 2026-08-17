"""Periodic repetition of a physical geometry."""
module Repeats

import StaticArrays: SVector
import ...Indices: ObstructionIndex
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: intersection_index_length, inside_index_length, contains_geometry_tuple
import ..Groups: inside_indices_for_any_type
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh, _translate_native
import ...Properties: GeometryProperties
import ..Groups

struct Repeat{N, P<:PhysicalGeometry{N}} <: PhysicalGeometry{N}
    geometry::P
    repeats::SVector{N, Float64}
    lower_overlap::SVector{N, Float64}
    upper_overlap::SVector{N, Float64}

    function Repeat{N, P}(geometry::P, repeats::SVector{N, Float64}) where {N, P<:PhysicalGeometry{N}}
        all(repeats .> 0) || throw(ArgumentError("repeat distances must be positive"))
        box = InternalBoundingBox(geometry)
        all(InternalBoundingBoxes.lower(box) .>= -1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far below its repeat bounds"))
        all(InternalBoundingBoxes.upper(box) .<= 1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far above its repeat bounds"))
        lower_overlap = max.(0.0, -repeats / 2 - InternalBoundingBoxes.lower(box))
        upper_overlap = max.(0.0, InternalBoundingBoxes.upper(box) - repeats / 2)
        new{N, P}(geometry, repeats, lower_overlap, upper_overlap)
    end
end

contains_geometry_tuple(::Type{<:Repeat{N, P}}) where {N, P} = contains_geometry_tuple(P)
intersection_index_length(::Type{<:Repeat{N, P}}) where {N, P} =
    contains_geometry_tuple(P) ? 0 : intersection_index_length(P)
inside_index_length(::Type{<:Repeat{N, P}}) where {N, P} =
    contains_geometry_tuple(P) ? 0 : inside_index_length(P)

Repeat(geometry::P, repeats::AbstractVector{<:Real}) where {N, P<:PhysicalGeometry{N}} =
    Repeat{N, P}(geometry, SVector{N, Float64}(repeats))

has_inside(::Type{<:Repeat{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Repeat}) = false

_wrap(repeat::Repeat, position) = mod.(position .+ repeat.repeats / 2, repeat.repeats) .- repeat.repeats / 2

_child_image(repeat::Repeat, position, image) = position .- image .* repeat.repeats

function _candidate_images(repeat::Repeat{N}, start, destination=start) where {N}
    local_start = _wrap(repeat, start)
    local_destination = _wrap(repeat, destination)
    child_box = InternalBoundingBox(repeat.geometry)
    child_lower = InternalBoundingBoxes.lower(child_box)
    child_upper = InternalBoundingBoxes.upper(child_box)
    images = SVector{N, Int}[]
    for image in Iterators.product(((-1):1 for _ in 1:N)...)
        image = SVector{N, Int}(image)
        child_start = _child_image(repeat, local_start, image)
        child_destination = _child_image(repeat, local_destination, image)
        lower = min.(child_start, child_destination)
        upper = max.(child_start, child_destination)
        all(lower .<= child_upper) && all(upper .>= child_lower) && push!(images, image)
    end
    images
end

function _needs_image_search(repeat::Repeat)
    child_box = InternalBoundingBox(repeat.geometry)
    child_lower = InternalBoundingBoxes.lower(child_box)
    child_upper = InternalBoundingBoxes.upper(child_box)
    any(
        ((repeat.lower_overlap .== 0) .& (repeat.upper_overlap .> 0) .&
         (child_upper .<= repeat.repeats)) .|
        ((repeat.lower_overlap .> 0) .& (repeat.upper_overlap .== 0) .&
         (child_lower .>= -repeat.repeats))
    )
end

function _candidate_shifts(repeat::Repeat{N}, start, destination=start) where {N}
    local_start = _wrap(repeat, start)
    local_destination = _wrap(repeat, destination)
    choices = ntuple(N) do dimension
        lower = repeat.lower_overlap[dimension] > 0 &&
            min(local_start[dimension], local_destination[dimension]) <=
                -repeat.repeats[dimension] / 2 + repeat.lower_overlap[dimension] &&
            max(local_start[dimension], local_destination[dimension]) >= -repeat.repeats[dimension] / 2
        upper = repeat.upper_overlap[dimension] > 0 &&
            min(local_start[dimension], local_destination[dimension]) <=
                repeat.repeats[dimension] / 2 &&
            max(local_start[dimension], local_destination[dimension]) >=
                repeat.repeats[dimension] / 2 - repeat.upper_overlap[dimension]
        values = [0]
        lower && push!(values, 1)
        upper && push!(values, -1)
        values
    end
    shifts = SVector{N, Int}[]
    for shift in Iterators.product(choices...)
        push!(shifts, SVector{N, Int}(shift))
    end
    shifts
end

function inside_indices(
    repeat::Repeat{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    local_position = _wrap(repeat, position)
    geometry_type = typeof(repeat.geometry)
    indices = contains_geometry_tuple(geometry_type) ?
        ObstructionIndex[] : ObstructionIndex{inside_index_length(geometry_type)}[]
    candidates = if Base.isempty(intersection) && _needs_image_search(repeat)
        ((_child_image(repeat, local_position, image), image) for image in _candidate_images(repeat, local_position))
    else
        ((local_position .+ shift .* repeat.repeats, shift) for shift in _candidate_shifts(repeat, local_position))
    end
    for (child_position, _) in candidates
        append!(indices, inside_indices_for_any_type(
            repeat.geometry,
            child_position,
            intersection,
        ))
    end
    sort!(unique!(indices), by=index -> Tuple(index.indices))
end

InternalBoundingBox(::Repeat) = throw(ArgumentError("repeated geometries do not have a finite bounding box"))

size_scale(repeat::Repeat) = min(size_scale(repeat.geometry), minimum(repeat.repeats))

function _repeat_samples(repeat::Repeat{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N}
    child_box = InternalBoundingBox(repeat.geometry)
    lower = floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats)
    upper = ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats)
    draws = (let
        displacement = SVector{N, Float64}(indices) .* repeat.repeats
        values = random_surface_positions(repeat.geometry, density,
            InternalBoundingBoxes.shift(bounding_box, -displacement), scale_density)
        ([value + displacement for value in values[1]], values[2])
    end for indices in Iterators.product((lower[index]:upper[index] for index in 1:N)...))
    Groups._combine(draws, Val(N))
end

random_surface_positions(repeat::Repeat, density::GeometryProperties, bounding_box, scale_density) =
    _repeat_samples(repeat, density, bounding_box, scale_density)

function _geometry_mesh(repeat::Repeat{N}; bounding_box=nothing, kwargs...) where N
    child_box = InternalBoundingBox(repeat.geometry)
    lower, upper = isnothing(bounding_box) ? (zeros(Int, N), zeros(Int, N)) :
        (floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats),
         ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats))
    child = _geometry_mesh(repeat.geometry; kwargs...)
    result = Any[]
    for index in Iterators.product((lower[i]:upper[i] for i in 1:N)...)
        displacement = SVector{N, Float64}(index) .* repeat.repeats
        if !isnothing(bounding_box)
            shifted_lower = InternalBoundingBoxes.lower(child_box) + displacement
            shifted_upper = InternalBoundingBoxes.upper(child_box) + displacement
            (all(shifted_lower .<= InternalBoundingBoxes.upper(bounding_box)) && all(shifted_upper .>= InternalBoundingBoxes.lower(bounding_box))) || continue
        end
        append!(result, [_translate_native(value, displacement) for value in child])
    end
    result
end

function detect_intersection(
    repeat::Repeat{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    displacement = destination - start
    iszero(displacement) && return Intersection{N, intersection_index_length(typeof(repeat))}()
    closest = Intersection{N, intersection_index_length(typeof(repeat))}()
    previous = previous_hit
    scaled_start = (start .+ repeat.repeats / 2) ./ repeat.repeats
    scaled_destination = (destination .+ repeat.repeats / 2) ./ repeat.repeats
    for (voxel, entry_time, _, exit_time, _) in ray_grid_intersections(scaled_start, scaled_destination)
        segment_start = start + entry_time .* displacement
        segment_destination = start + exit_time .* displacement
        cell_shift = voxel .* repeat.repeats
        local_start = segment_start - cell_shift
        local_destination = segment_destination - cell_shift
        for shift in _candidate_shifts(repeat, local_start, local_destination)
            local_intersection = detect_intersection(
                repeat.geometry,
                local_start .+ shift .* repeat.repeats,
                local_destination .+ shift .* repeat.repeats,
                previous,
            )
            previous = Intersection{3, intersection_index_length(typeof(repeat))}()
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
        if !isempty(closest)
            return closest
        end
    end
    return Intersection{N, intersection_index_length(typeof(repeat))}()
end

end
