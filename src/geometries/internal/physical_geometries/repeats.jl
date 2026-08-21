"""Periodic repetition of a physical geometry."""
module Repeats

import StaticArrays: SVector
import Random: rand
import ...InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, get_child, has_inside, has_single_inside, inside_indices_eltype, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, distance_to_surface, _geometry_mesh, _translate_native, to_property_index
import ...Properties: GeometryProperties
import ..Groups
import ..Transformations: Shift

struct Repeat{N, P<:PhysicalGeometry{N}} <: Groups.GroupGeometry{N, Shift{N, P}}
    geometry::P
    repeats::SVector{N, Float64}
    normalized_bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}

    function Repeat{N, P}(geometry::P, repeats::SVector{N, Float64}) where {N, P<:PhysicalGeometry{N}}
        all(repeats .> 0) || throw(ArgumentError("repeat distances must be positive"))
        box = InternalBoundingBox(geometry)
        all(InternalBoundingBoxes.lower(box) .>= -1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far below its repeat bounds"))
        all(InternalBoundingBoxes.upper(box) .<= 1.5 .* repeats) ||
            throw(ArgumentError("geometry extends too far above its repeat bounds"))
        normalized_lower = InternalBoundingBoxes.lower(box) ./ repeats
        normalized_upper = InternalBoundingBoxes.upper(box) ./ repeats
        normalized_bounding_box = InternalBoundingBoxes.InternalBoundingBox{N}(
            (normalized_upper - normalized_lower) / 2,
            (normalized_upper + normalized_lower) / 2,
        )
        new{N, P}(geometry, repeats, normalized_bounding_box)
    end
end

function Base.show(io::IO, ::Type{T}) where {N, P, T <: Repeat{N, P}}
    print(io, "Repeat{")
    show(io, P)
    print(io, "}")
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
    iszero(destination - start) && return ()
    direction = SVector{N, Bool}(destination .>= start)
    start_repeat = directed_repeat(repeat, start, direction)
    destination_repeat = directed_repeat(repeat, destination, .!direction)
    ranges = (
        min(start_repeat[i], destination_repeat[i]):max(start_repeat[i], destination_repeat[i])
        for i in 1:N
    )
    (
        let copy_shift = SVector{N, Int}(repeat_index)
            (copy_shift, Shift(repeat.geometry, copy_shift .* repeat.repeats), 0.0)
        end
        for repeat_index in Iterators.product(ranges...)
    )
end

function find_intersection(
    repeat::Repeat{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
) where {N}
    current = nothing
    direction = SVector{N, Bool}(destination .>= start)
    start_repeat = directed_repeat(repeat, start, direction)
    destination_repeat = directed_repeat(repeat, destination, .!direction)
    ranges = (
        min(start_repeat[i], destination_repeat[i]):max(start_repeat[i], destination_repeat[i])
        for i in 1:N
    )
    for repeat_index in Iterators.product(ranges...)
        copy_shift = SVector{N, Int}(repeat_index)
        displacement = copy_shift .* repeat.repeats
        candidate_previous_hit = isnothing(previous_hit) || previous_hit[1] != copy_shift ?
            nothing : previous_hit[2:end]
        intersect = find_intersection(
            repeat.geometry,
            start - displacement,
            destination - displacement,
            candidate_previous_hit,
        )
        isnothing(intersect) && continue
        candidate = (copy_shift, intersect...)
        if isnothing(current) || candidate[end] < current[end]
            current = candidate
        end
    end
    current
end

function get_child(repeat::Repeat{N}, indices::Tuple) where {N}
    copy_shift = indices[1]
    Shift(repeat.geometry, copy_shift .* repeat.repeats), indices[2:end]
end

function Groups.inside_candidates(
    repeat::Repeat{N},
    position::SVector{N, Float64},
) where {N}
    (
        let copy_shift = SVector{N, Int}(repeat_index)
            (copy_shift, Shift(repeat.geometry, copy_shift .* repeat.repeats))
        end
        for repeat_index in Iterators.product(
            (first_repeat(repeat, position)[i]:last_repeat(repeat, position)[i] for i in 1:N)...
        )
    )
end

function random_surface_positions(
    repeat::Repeat{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N},
    scale_density,
) where {N}
    lower = first_repeat(repeat, InternalBoundingBoxes.lower(bounding_box))
    upper = last_repeat(repeat, InternalBoundingBoxes.upper(bounding_box))
    any(lower .> upper) && return SVector{N, Float64}[], Tuple[]

    ranges = ntuple(index -> lower[index]:upper[index], N)
    nrepeats = prod(length, ranges)
    positions, indices = random_surface_positions(
        repeat.geometry,
        density,
        InternalBoundingBox(repeat.geometry),
        scale_density * nrepeats,
    )
    repeat_indices = (
        SVector{N, Int}(rand(ranges[index]) for index in 1:N)
        for _ in eachindex(positions)
    )
    shifted_positions = [
        position + repeat_index .* repeat.repeats
        for (position, repeat_index) in zip(positions, repeat_indices)
    ]
    shifted_indices = [
        (repeat_index, child_indices...)
        for (repeat_index, child_indices) in zip(repeat_indices, indices)
    ]
    shifted_positions, shifted_indices
end

has_inside(::Type{<:Repeat{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Repeat}) = false

function first_repeat(repeat::Repeat{N}, position::SVector{N, Float64}) where {N}
    Int.(ceil.((position ./ repeat.repeats) .-
        InternalBoundingBoxes.upper(repeat.normalized_bounding_box)))
end

function last_repeat(repeat::Repeat{N}, position::SVector{N, Float64}) where {N}
    Int.(floor.((position ./ repeat.repeats) .-
        InternalBoundingBoxes.lower(repeat.normalized_bounding_box)))
end

function directed_repeat(
    repeat::Repeat{N},
    position::SVector{N, Float64},
    direction::SVector{N, Bool},
) where {N}
    first = first_repeat(repeat, position)
    last = last_repeat(repeat, position)
    SVector{N, Int}(direction[i] ? first[i] : last[i] for i in 1:N)
end

function _repeat_indices(::Repeat{N}, lower, upper) where {N}
    Iterators.product((lower[i]:upper[i] for i in 1:N)...)
end

function distance_to_surface(
    repeat::Repeat{N},
    position::SVector{N, Float64},
) where {N}
    local_position = position .- repeat.repeats .* round.(position ./ repeat.repeats)
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

InternalBoundingBox(::Repeat) = throw(ArgumentError("repeated geometries do not have a finite bounding box"))

size_scale(repeat::Repeat) = min(size_scale(repeat.geometry), minimum(repeat.repeats))

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
