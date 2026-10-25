"""Grid acceleration over the leaf surfaces of a finite geometry tree."""
module Grids

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, InternalBoundingBox, leaf_bounding_boxes
import ..PhysicalGeometries: intersection_type, inside_indices_eltype, find_intersection
import ..PhysicalGeometries: rewrap, has_equivalent_wrapper
import ..PhysicalGeometries: _volume_inside_indices, to_inside_index
import ..Transparents: Transparent
import ..GridDispatch: IntersectionGrid, GridIterator
import ...InternalBoundingBoxes

"""
    Grid(geometry::PhysicalGeometry; resolution=nothing)

Wrap a finite geometry in an intersection grid over its leaf bounding boxes.
`resolution` specifies the cell size; `nothing` uses the `IntersectionGrid`
estimate. The requested value is stored in `user_resolution` and reused when
rewrapping a simplified child. Leaf indices follow the hierarchy returned by
`leaf_bounding_boxes`, so intersection paths and properties retain their original
geometry indices.

The wrapper expects an unrestricted incoming candidate and supplies individual
candidate paths to its child during grid traversal. For periodic geometries,
place the grid inside `Repeat`.
"""
struct Grid{N, P<:PhysicalGeometry{N}, I, V} <: Transparent{N, P}
    geometry::P
    grid::IntersectionGrid{N}
    indices::Vector{I}
    fully_inside_indices::Array{Vector{V}, N}
    user_resolution::Union{Nothing, Float64}
end

function Grid(geometry::P; resolution=nothing) where {N, P<:PhysicalGeometry{N}}
    boxes = InternalBoundingBox{N}[]
    indices = Vector{intersection_type(P)}()
    _flatten_leaf_bounding_boxes!(boxes, indices, leaf_bounding_boxes(geometry), ())
    grid = IntersectionGrid(boxes; resolution)
    inside_indices_by_voxel = _fully_inside_indices(geometry, grid, indices)
    Grid{N, P, eltype(indices), inside_indices_eltype(P)}(
        geometry, grid, indices, inside_indices_by_voxel, resolution,
    )
end

function _fully_inside_indices(geometry, grid::IntersectionGrid{N}, collision_indices) where {N}
    inside_type = inside_indices_eltype(typeof(geometry))
    result = Array{Vector{inside_type}, N}(undef, size(grid.indices))
    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    cell_size = 1 ./ grid.inv_resolution
    for voxel in CartesianIndices(result)
        centre = lower_bound .+ (SVector{N, Float64}(Tuple(voxel)) .- 0.5) .* cell_size
        inside_at_centre = _volume_inside_indices(geometry, centre)
        collision_inside = Set{inside_type}()
        for collision_index in grid.indices[voxel]
            inside_index = to_inside_index(geometry, collision_indices[collision_index])
            isnothing(inside_index) || push!(collision_inside, inside_index)
        end
        result[voxel] = inside_type[
            inside_index for inside_index in inside_at_centre
            if !(inside_index in collision_inside)
        ]
    end
    result
end

rewrap(wrapper::Grid{N}, geometry::PhysicalGeometry{N}) where {N} =
    Grid(geometry; resolution=wrapper.user_resolution)
has_equivalent_wrapper(first::Grid{N}, second::Grid{N}) where {N} =
    first.user_resolution == second.user_resolution

function _flatten_leaf_bounding_boxes!(boxes, indices, box::InternalBoundingBox, path)
    push!(boxes, box)
    push!(indices, path)
    nothing
end

function _flatten_leaf_bounding_boxes!(boxes, indices, children::Union{Tuple, AbstractVector}, path)
    for (index, child) in enumerate(children)
        _flatten_leaf_bounding_boxes!(boxes, indices, child, (path..., index))
    end
    nothing
end

InternalBoundingBox(geometry::Grid; kwargs...) = geometry.grid.bounding_box

function _earlier_hit(hit, best)
    isnothing(best) && return true
    hit[end] < best[end] && return true
    hit[end] == best[end] && isless(hit[1:(end - 2)], best[1:(end - 2)])
end

function find_intersection(
    geometry::Grid{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
    inside=nothing,
    candidate::Union{Nothing, Tuple}=nothing,
) where {N}
    isnothing(candidate) || throw(ArgumentError("Grid requires an unrestricted incoming candidate"))
    best = nothing
    for (index, entry_time) in GridIterator(geometry.grid, start, destination)
        !isnothing(best) && best[end] < entry_time && break
        hit = find_intersection(
            geometry.geometry, start, destination, previous_hit, inside, geometry.indices[index],
        )
        if !isnothing(hit) && _earlier_hit(hit, best)
            best = hit
        end
    end
    best
end

end
