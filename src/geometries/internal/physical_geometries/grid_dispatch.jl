module GridDispatch

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: finalise!

struct IntersectionGrid{N}
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    grid_bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    inv_resolution::SVector{N, Float64}
    dimensions::SVector{N, Int}
    number_of_children::Int
    indices::Array{Vector{Int}, N}
end

"""Iterate over grid child indices in ray order."""
struct GridIterator{N, I, R}
    indices::I
    ray::R
    grid_size::SVector{N, Int}
end

function GridIterator(
    grid::IntersectionGrid{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
) where {N}
    indices = _require_grid_indices(grid)
    grid_size = grid.dimensions
    iszero(destination - start) && return GridIterator{N, typeof(indices), Nothing}(indices, nothing, grid_size)
    InternalBoundingBoxes.does_intersect(grid.grid_bounding_box, start, destination) ||
        return GridIterator{N, typeof(indices), Nothing}(indices, nothing, grid_size)

    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    scaled_start = (start - lower_bound) .* grid.inv_resolution
    scaled_destination = (destination - lower_bound) .* grid.inv_resolution
    GridIterator{N, typeof(indices), typeof(ray_grid_intersections(scaled_start, scaled_destination))}(
        indices,
        ray_grid_intersections(scaled_start, scaled_destination),
        grid_size,
    )
end

Base.IteratorSize(::Type{<:GridIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:GridIterator}) = Tuple{Int, Float64}

Base.iterate(iterator::GridIterator{N, I, Nothing}) where {N, I} = nothing

function Base.iterate(iterator::GridIterator)
    state = _first_grid_state(iterator)
    isnothing(state) && return nothing
    _iterate_grid_candidates(iterator, state)
end

Base.iterate(iterator::GridIterator, state) = _iterate_grid_candidates(iterator, state)

function _first_grid_state(iterator::GridIterator)
    next_voxel = iterate(iterator.ray)
    while !isnothing(next_voxel)
        ((voxel, entry_time, _, _), ray_state) = next_voxel
        if !any(voxel .< 0) && !any(voxel .>= iterator.grid_size)
            children = iterator.indices[voxel .+ 1...]
            return (ray_state, children, 1, true, entry_time)
        end
        next_voxel = iterate(iterator.ray, ray_state)
    end
    nothing
end

function _iterate_grid_candidates(iterator::GridIterator, state)
    ray_state, children, child_position, entered_grid, entry_time = state

    while true
        while child_position <= length(children)
            child_index = children[child_position]
            child_position += 1
            return (
                (child_index, entry_time),
                (ray_state, children, child_position, entered_grid, entry_time),
            )
        end

        next_voxel = isnothing(ray_state) ? iterate(iterator.ray) : iterate(iterator.ray, ray_state)
        isnothing(next_voxel) && return nothing
        ((voxel, entry_time, _, _), ray_state) = next_voxel

        if any(voxel .< 0) || any(voxel .>= iterator.grid_size)
            entered_grid && return nothing
            continue
        end

        entered_grid = true
        children = iterator.indices[voxel .+ 1...]
        child_position = 1
    end
end

function _require_grid_indices(grid::IntersectionGrid)
    isassigned(grid.indices, lastindex(grid.indices)) || throw(ArgumentError(
        "IntersectionGrid indices are delayed; call finalise!(grid, bounding_boxes) before traversal",
    ))
    grid.indices
end

function _grid_bounding_box(box::InternalBoundingBoxes.InternalBoundingBox{N}, minimum_thickness) where {N}
    lower_bound = InternalBoundingBoxes.lower(box)
    upper_bound = InternalBoundingBoxes.upper(box)
    extent = upper_bound - lower_bound
    grid_extent = max.(extent, minimum_thickness)
    InternalBoundingBoxes.InternalBoundingBox(
        grid_extent / 2,
        (lower_bound + upper_bound) / 2,
    )
end

"""Build an intersection grid; `delayed=true` defers its per-cell indices."""
function IntersectionGrid(
    bounding_boxes::AbstractVector{<:InternalBoundingBoxes.InternalBoundingBox{N}};
    resolution=nothing,
    bounding_box=nothing,
    delayed::Bool=false,
) where {N}
    isempty(bounding_boxes) && throw(ArgumentError("cannot construct a grid for an empty bounding-box vector"))
    bounding_box = if isnothing(bounding_box)
        lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in bounding_boxes))
        upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in bounding_boxes))
        InternalBoundingBoxes.InternalBoundingBox(
            (upper_bound - lower_bound) / 2,
            (upper_bound + lower_bound) / 2,
        )
    else
        bounding_box
    end
    extent = InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(bounding_box)
    maximum_extent = maximum(extent)
    resolution = isnothing(resolution) ?
        (iszero(maximum_extent) ? 1.0 : maximum_extent / max(ceil(Int, length(bounding_boxes)^(1 / N)), 1)) : resolution
    isinf(resolution) || resolution > 0 || throw(ArgumentError("grid resolution must be positive"))
    # Keep single-cell (Inf) grids finite, including fully zero-extent geometry.
    minimum_thickness = isinf(resolution) ? (iszero(maximum_extent) ? 1.0 : maximum_extent) : resolution
    grid_bounding_box = _grid_bounding_box(bounding_box, minimum_thickness)
    extent = InternalBoundingBoxes.upper(grid_bounding_box) - InternalBoundingBoxes.lower(grid_bounding_box)
    dimensions = isinf(resolution) ? SVector{N, Int}(ones(Int, N)) :
        max.(SVector{N, Int}(ceil.(Int, extent ./ resolution)), 1)
    cell_size = extent ./ dimensions
    indices = delayed ? Array{Vector{Int}, N}(undef, Tuple(dimensions)) :
        InternalBoundingBoxes.grid_indices(grid_bounding_box, dimensions, bounding_boxes)
    IntersectionGrid{N}(
        bounding_box,
        grid_bounding_box,
        SVector{N, Float64}(1 ./ cell_size),
        dimensions,
        length(bounding_boxes),
        indices,
    )
end

"""Populate the per-cell candidate indices of a delayed intersection grid."""
function finalise!(
    grid::IntersectionGrid{N},
    bounding_boxes::AbstractVector{<:InternalBoundingBoxes.InternalBoundingBox{N}},
) where {N}
    length(bounding_boxes) == grid.number_of_children || throw(ArgumentError(
        "cannot finalise IntersectionGrid with a different number of bounding boxes",
    ))
    isassigned(grid.indices, lastindex(grid.indices)) && return grid
    indices = InternalBoundingBoxes.grid_indices(
        grid.grid_bounding_box, grid.dimensions, bounding_boxes,
    )
    for index in eachindex(indices)
        grid.indices[index] = indices[index]
    end
    grid
end

end
