module GridDispatch

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: PhysicalGeometry, filter_bounding_box, finalise!

struct IntersectionGrid{N}
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    grid_bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    inv_resolution::SVector{N, Float64}
    dimensions::SVector{N, Int}
    number_of_children::Int
    geometries::Array{Union{Nothing, PhysicalGeometry{N}}, N}
end

"""Iterate over non-empty cell geometries in ray order."""
struct GridIterator{N, G, R}
    geometries::G
    ray::R
    grid_size::SVector{N, Int}
end

function GridIterator(
    grid::IntersectionGrid{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
) where {N}
    geometries = _require_grid_geometries(grid)
    grid_size = grid.dimensions
    iszero(destination - start) && return GridIterator{N, typeof(geometries), Nothing}(geometries, nothing, grid_size)
    InternalBoundingBoxes.does_intersect(grid.grid_bounding_box, start, destination) ||
        return GridIterator{N, typeof(geometries), Nothing}(geometries, nothing, grid_size)

    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    scaled_start = (start - lower_bound) .* grid.inv_resolution
    scaled_destination = (destination - lower_bound) .* grid.inv_resolution
    GridIterator{N, typeof(geometries), typeof(ray_grid_intersections(scaled_start, scaled_destination))}(
        geometries,
        ray_grid_intersections(scaled_start, scaled_destination),
        grid_size,
    )
end

Base.IteratorSize(::Type{<:GridIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:GridIterator{N}}) where {N} = Tuple{PhysicalGeometry{N}, Float64}

Base.iterate(iterator::GridIterator{N, G, Nothing}) where {N, G} = nothing

function Base.iterate(iterator::GridIterator)
    _iterate_grid_geometries(iterator, nothing, false)
end

Base.iterate(iterator::GridIterator, state) = _iterate_grid_geometries(iterator, state...)

function _iterate_grid_geometries(iterator::GridIterator, ray_state, entered_grid)
    next_voxel = isnothing(ray_state) ? iterate(iterator.ray) : iterate(iterator.ray, ray_state)
    while !isnothing(next_voxel)
        ((voxel, entry_time, _, _), ray_state) = next_voxel
        if any(voxel .< 0) || any(voxel .>= iterator.grid_size)
            entered_grid && return nothing
        else
            entered_grid = true
            geometry = iterator.geometries[voxel .+ 1...]
            if !isnothing(geometry)
                return ((geometry, entry_time), (ray_state, entered_grid))
            end
        end
        next_voxel = iterate(iterator.ray, ray_state)
    end
    nothing
end

function _require_grid_geometries(grid::IntersectionGrid)
    isassigned(grid.geometries, lastindex(grid.geometries)) || throw(ArgumentError(
        "IntersectionGrid is not finalised; call finalise!(grid, geometry) before traversal",
    ))
    grid.geometries
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

function cell_bounding_box(grid::IntersectionGrid{N}, coordinate::CartesianIndex{N}) where {N}
    cell_size = inv.(grid.inv_resolution)
    lower = InternalBoundingBoxes.lower(grid.grid_bounding_box) .+
        (SVector{N, Float64}(Tuple(coordinate)) .- 1) .* cell_size
    upper = lower + cell_size
    padding = 30 * eps(max(1.0, maximum(abs.(lower)), maximum(abs.(upper))))
    InternalBoundingBoxes.InternalBoundingBox(
        cell_size / 2 .+ padding,
        lower + cell_size / 2,
    )
end

"""Build an unfinalised spatial grid layout from leaf bounding boxes."""
function IntersectionGrid(
    bounding_boxes::AbstractVector{<:InternalBoundingBoxes.InternalBoundingBox{N}};
    resolution=nothing,
    bounding_box=nothing,
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
    geometries = Array{Union{Nothing, PhysicalGeometry{N}}, N}(undef, Tuple(dimensions))
    IntersectionGrid{N}(
        bounding_box,
        grid_bounding_box,
        SVector{N, Float64}(1 ./ cell_size),
        dimensions,
        length(bounding_boxes),
        geometries,
    )
end

"""Populate the per-cell geometry payloads by filtering `geometry` to each cell."""
function finalise!(
    grid::IntersectionGrid{N}, geometry::PhysicalGeometry{N},
) where {N}
    isassigned(grid.geometries, lastindex(grid.geometries)) && return grid
    for coordinate in CartesianIndices(grid.geometries)
        grid.geometries[coordinate] = filter_bounding_box(
            geometry, cell_bounding_box(grid, coordinate),
        )
    end
    grid
end

end
