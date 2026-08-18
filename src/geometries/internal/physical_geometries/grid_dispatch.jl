module GridDispatch

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: Intersection

struct IntersectionGrid{N}
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    grid_bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    inv_resolution::SVector{N, Float64}
    indices::Array{Vector{Int}, N}
end

function _grid_bounding_box(box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    lower_bound = InternalBoundingBoxes.lower(box)
    upper_bound = InternalBoundingBoxes.upper(box)
    extent = upper_bound - lower_bound
    grid_extent = max.(extent, eps(Float64))
    InternalBoundingBoxes.InternalBoundingBox(
        grid_extent / 2,
        (lower_bound + upper_bound) / 2,
    )
end

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
    grid_bounding_box = _grid_bounding_box(bounding_box)
    extent = InternalBoundingBoxes.upper(grid_bounding_box) - InternalBoundingBoxes.lower(grid_bounding_box)
    resolution = isnothing(resolution) ?
        maximum(extent) / max(ceil(Int, length(bounding_boxes)^(1 / N)), 1) : resolution
    isinf(resolution) || resolution > 0 || throw(ArgumentError("grid resolution must be positive"))
    dimensions = isinf(resolution) ? SVector{N, Int}(ones(Int, N)) :
        max.(SVector{N, Int}(ceil.(Int, extent ./ resolution)), 1)
    cell_size = extent ./ dimensions
    IntersectionGrid{N}(
        bounding_box,
        grid_bounding_box,
        SVector{N, Float64}(1 ./ cell_size),
        InternalBoundingBoxes.grid_indices(grid_bounding_box, dimensions, bounding_boxes),
    )
end

"""Dispatch a ray through a precomputed grid and stop at the first hit voxel."""
function detect_intersection_grid(
    detect_child,
    grid::IntersectionGrid{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3},
    empty_intersection::Intersection{N, M},
) where {N, M}
    iszero(destination - start) && return empty_intersection
    InternalBoundingBoxes.could_intersect(grid.grid_bounding_box, start, destination) || return empty_intersection
    InternalBoundingBoxes.does_intersect(grid.grid_bounding_box, start, destination) || return empty_intersection

    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    scaled_start = (start - lower_bound) .* grid.inv_resolution
    scaled_destination = (destination - lower_bound) .* grid.inv_resolution
    found = empty_intersection
    entered_grid = false
    grid_size = SVector{N, Int}(size(grid.indices))

    for (voxel, _, _, exit_time, _) in ray_grid_intersections(scaled_start, scaled_destination)
        if any(voxel .< 0) || any(voxel .>= grid_size)
            entered_grid && break
            continue
        end
        entered_grid = true

        for child_index in grid.indices[voxel .+ 1...]
            intersection = detect_child(child_index, previous_hit)
            if intersection.distance < found.distance
                found = intersection
            end
        end
        found.distance <= exit_time && return found
    end
    found
end

end
