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
    number_of_leaves::Int
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

function _range_bounding_box(grid::IntersectionGrid{N}, ranges::NTuple{N, UnitRange{Int}}) where {N}
    cell_size = inv.(grid.inv_resolution)
    grid_lower = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    lower_coordinate = SVector{N, Float64}(ntuple(d -> first(ranges[d]) - 1, N))
    upper_coordinate = SVector{N, Float64}(ntuple(d -> last(ranges[d]), N))
    lower = grid_lower + lower_coordinate .* cell_size
    upper = grid_lower + upper_coordinate .* cell_size
    padding = 30 * eps(max(1.0, maximum(abs.(lower)), maximum(abs.(upper))))
    InternalBoundingBoxes.InternalBoundingBox(
        (upper - lower) / 2 .+ padding,
        (lower + upper) / 2,
    )
end

function cell_bounding_box(grid::IntersectionGrid{N}, coordinate::CartesianIndex{N}) where {N}
    _range_bounding_box(grid, ntuple(d -> coordinate[d]:coordinate[d], N))
end

function _populate_grid_range!(grid::IntersectionGrid{N}, geometry, ranges) where {N}
    if isnothing(geometry)
        for coordinate in CartesianIndices(ranges)
            grid.geometries[coordinate] = nothing
        end
        return grid
    end

    split_dimension = argmax(length.(ranges))
    if length(ranges[split_dimension]) == 1
        coordinate = CartesianIndex(ntuple(d -> first(ranges[d]), N))
        grid.geometries[coordinate] = filter_bounding_box(
            geometry, _range_bounding_box(grid, ranges),
        )
        return grid
    end

    range = ranges[split_dimension]
    midpoint = first(range) + fld(length(range) - 1, 2)
    for child_range in (first(range):midpoint, (midpoint + 1):last(range))
        child_ranges = ntuple(
            d -> d == split_dimension ? child_range : ranges[d], N,
        )
        child_geometry = filter_bounding_box(
            geometry, _range_bounding_box(grid, child_ranges),
        )
        _populate_grid_range!(grid, child_geometry, child_ranges)
    end
    grid
end

"""Build an unfinalised grid layout from the full bounds and leaf count."""
function IntersectionGrid(
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N},
    number_of_leaves::Integer;
    resolution=nothing,
) where {N}
    number_of_leaves > 0 || throw(ArgumentError("cannot construct a grid for geometry without leaves"))
    extent = InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(bounding_box)
    maximum_extent = maximum(extent)
    resolution = isnothing(resolution) ?
        (iszero(maximum_extent) ? 1.0 : maximum_extent / max(ceil(Int, number_of_leaves^(1 / N)), 1)) : resolution
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
        Int(number_of_leaves),
        geometries,
    )
end

"""Populate the per-cell geometry payloads by filtering `geometry` to each cell."""
function finalise!(
    grid::IntersectionGrid{N}, geometry::PhysicalGeometry{N},
) where {N}
    isassigned(grid.geometries, lastindex(grid.geometries)) && return grid
    ranges = ntuple(d -> 1:grid.dimensions[d], N)
    _populate_grid_range!(grid, geometry, ranges)
end

end
