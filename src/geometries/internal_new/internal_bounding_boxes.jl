"""Compact internal representations of axis-aligned bounding boxes."""
module InternalBoundingBoxes

import StaticArrays: SVector


"""
    InternalBoundingBox{N}(half_size, center=Nothing)

Bounding box used internally by the geometry engine. Bounds are normalized to
`lower` and `upper` at construction time for efficient intersection testing.
"""
struct InternalBoundingBox{N}
    lower::SVector{N, Float64}
    upper::SVector{N, Float64}

    function InternalBoundingBox{N}(half_size::SVector{N, Float64}, center::SVector{N, Float64}) where {N}
        all(half_size .>= 0) || throw(ArgumentError("bounding-box half-sizes must be nonnegative"))
        new{N}(center - half_size, center + half_size)
    end
end


function InternalBoundingBox{N}(
    half_size::Union{Real, AbstractVector},
    center::Union{Nothing, AbstractVector}=nothing,
) where {N}
    half_size_int = half_size isa Real ? SVector{N, Float64}(fill(half_size, N)) : SVector{N, Float64}(half_size)
    center_int = isnothing(center) ? zero(SVector{N, Float64}) : SVector{N, Float64}(center)
    return InternalBoundingBox{N}(half_size_int, center_int)
end

InternalBoundingBox(half_size::AbstractVector) = InternalBoundingBox{length(half_size)}(half_size)
InternalBoundingBox(half_size, center::AbstractVector) = InternalBoundingBox{length(center)}(half_size, center)


lower(box::InternalBoundingBox) = box.lower
upper(box::InternalBoundingBox) = box.upper
center(box::InternalBoundingBox) = (box.lower + box.upper) / 2
half_size(box::InternalBoundingBox) = (box.upper - box.lower) / 2


"""
    grid_indices(box, dimensions, bounding_boxes)

Return the indices of the grid cells intersected by each bounding box. The grid
does not contain periodic copies of the bounding boxes.
"""
function grid_indices(
    box::InternalBoundingBox{N},
    dimensions::AbstractVector{<:Integer},
    bounding_boxes::AbstractVector{<:InternalBoundingBox{N}},
) where N
    dimensions = SVector{N, Int}(dimensions)
    all(dimensions .> 0) || throw(ArgumentError("grid dimensions must be positive"))

    cell_size = (upper(box) - lower(box)) ./ dimensions
    indices = Array{Vector{Int}, N}(undef, Tuple(dimensions))
    for coordinate in CartesianIndices(indices)
        indices[coordinate] = Int[]
    end

    for (box_index, child_box) in enumerate(bounding_boxes)
        lower_coordinate = max.(
            Int.(floor.((lower(child_box) - lower(box)) ./ cell_size)) .+ 1,
            1,
        )
        upper_coordinate = min.(
            Int.(ceil.((upper(child_box) - lower(box)) ./ cell_size)),
            dimensions,
        )
        all(lower_coordinate .<= upper_coordinate) || continue
        for coordinate in Iterators.product(
            (lower_coordinate[i]:upper_coordinate[i] for i in 1:N)...,
        )
            push!(indices[coordinate...], box_index)
        end
    end
    indices
end


shift(box::InternalBoundingBox, displacement::AbstractVector{<:Real}) = begin
    displacement_int = SVector{length(box.lower), Float64}(displacement)
    InternalBoundingBox{length(box.lower)}(half_size(box), center(box) + displacement_int)
end

isinside(box::InternalBoundingBox, position::AbstractVector) =
    all(position .>= box.lower) && all(position .<= box.upper)

function could_intersect(box::InternalBoundingBox{N}, start::AbstractVector, destination::AbstractVector) where {N}
    return all(
        (start[i] <= box.upper[i] || destination[i] <= box.upper[i]) &&
        (start[i] >= box.lower[i] || destination[i] >= box.lower[i])
        for i in 1:N
    )
end

function does_intersect(box::InternalBoundingBox{N}, start::AbstractVector, destination::AbstractVector) where {N}
    t_min = 0.0
    t_max = 1.0

    for i in 1:N
        displacement = destination[i] - start[i]
        if iszero(displacement)
            box.lower[i] <= start[i] <= box.upper[i] || return false
        else
            t1 = (box.lower[i] - start[i]) / displacement
            t2 = (box.upper[i] - start[i]) / displacement
            t_min = max(t_min, min(t1, t2))
            t_max = min(t_max, max(t1, t2))
            t_min <= t_max || return false
        end
    end
    return true
end

end
