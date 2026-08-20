"""Periodic repetition of a physical geometry."""
module Repeats

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, get_child, has_inside, has_single_inside, inside_indices_eltype, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, distance_to_surface, _geometry_mesh, _translate_native, to_property_index
import ...Properties: GeometryProperties
import ..Groups
import ..Transformations: Shift

struct Repeat{N, P<:PhysicalGeometry{N}} <: Groups.GroupGeometry{N, Shift{N, P}}
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

inside_indices_eltype(::Type{<:Repeat{N, P}}) where {N, P} =
    Groups._prepend_type(SVector{N, Int}, inside_indices_eltype(child_type(Repeat{N, P})))

Repeat(geometry::P, repeats::AbstractVector{<:Real}) where {N, P<:PhysicalGeometry{N}} =
    Repeat{N, P}(geometry, SVector{N, Float64}(repeats))

function to_property_index(geometry::Repeat, indices)
    if isnothing(indices)
        return nothing
    end
    child, child_indices = get_child(geometry, indices)
    cleaned = to_property_index(child, child_indices)
    return cleaned
end

function Groups.intersection_candidates(
    repeat::Repeat{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
) where {N}
    displacement = destination - start
    iszero(displacement) && return ()

    scaled_start = (start .+ repeat.repeats / 2) ./ repeat.repeats
    scaled_destination = (destination .+ repeat.repeats / 2) ./ repeat.repeats
    checked = Set{SVector{N, Int}}()

    (
        let copy_shift = voxel - local_shift
            (copy_shift, Shift(repeat.geometry, copy_shift .* repeat.repeats), entry_time)
        end
        for (voxel, entry_time, _, exit_time, _) in ray_grid_intersections(scaled_start, scaled_destination)
        for local_shift in _candidate_shifts(
            repeat,
            start + entry_time .* displacement - voxel .* repeat.repeats,
            start + exit_time .* displacement - voxel .* repeat.repeats,
        )
        if !((voxel - local_shift) in checked) && (push!(checked, voxel - local_shift); true)
    )
end

function get_child(repeat::Repeat{N}, indices::Tuple) where {N}
    copy_shift = indices[1]
    Shift(repeat.geometry, copy_shift .* repeat.repeats), indices[2:end]
end

function Groups.inside_candidates(
    repeat::Repeat{N},
    position::SVector{N, Float64},
) where {N}
    local_position = _wrap(repeat, position)
    cell_shift = SVector{N, Int}(round.(Int, (position - local_position) ./ repeat.repeats))
    checked = Set{SVector{N, Int}}()

    (
        let copy_shift = cell_shift - local_shift
            (copy_shift, Shift(repeat.geometry, copy_shift .* repeat.repeats))
        end
        for local_shift in _candidate_shifts(repeat, local_position)
        if !((cell_shift - local_shift) in checked) && (push!(checked, cell_shift - local_shift); true)
    )
end

has_inside(::Type{<:Repeat{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Repeat}) = false

_wrap(repeat::Repeat, position) = mod.(position .+ repeat.repeats / 2, repeat.repeats) .- repeat.repeats / 2

function distance_to_surface(
    repeat::Repeat{N},
    position::SVector{N, Float64},
) where {N}
    local_position = _wrap(repeat, position)
    minimum(
        (
            distance_to_surface(
                repeat.geometry,
                local_position - SVector{N, Float64}(shift) .* repeat.repeats,
            )
            for shift in Iterators.product(ntuple(_ -> -1:1, N)...)
        );
        init=Inf,
    )
end

function _candidate_shifts(repeat::Repeat{N}, start, destination=start) where {N}
    local_start = _wrap(repeat, start)
    local_destination = _wrap(repeat, destination)
    box = InternalBoundingBox(repeat.geometry)
    choices = ntuple(N) do dimension
        local_lower = min(local_start[dimension], local_destination[dimension])
        local_upper = max(local_start[dimension], local_destination[dimension])
        first_shift = ceil(Int, (InternalBoundingBoxes.lower(box)[dimension] - local_upper) / repeat.repeats[dimension])
        last_shift = floor(Int, (InternalBoundingBoxes.upper(box)[dimension] - local_lower) / repeat.repeats[dimension])
        shifts = collect(first_shift:last_shift)
        [0, filter(!iszero, shifts)...]
    end
    shifts = SVector{N, Int}[]
    for shift in Iterators.product(choices...)
        push!(shifts, SVector{N, Int}(shift))
    end
    shifts
end

InternalBoundingBox(::Repeat) = throw(ArgumentError("repeated geometries do not have a finite bounding box"))

size_scale(repeat::Repeat) = min(size_scale(repeat.geometry), minimum(repeat.repeats))

function Groups.group_geometries_bounded(repeat::Repeat{N}, bounding_box::InternalBoundingBox{N}; kwargs...) where {N}
    child_box = InternalBoundingBox(repeat.geometry)
    lower = floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats)
    upper = ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats)
    (
        let copy_shift = SVector{N, Int}(shift)
            (copy_shift, Shift(repeat.geometry, copy_shift .* repeat.repeats))
        end
        for shift in Iterators.product((lower[index]:upper[index] for index in 1:N)...)
    )
end

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

end
