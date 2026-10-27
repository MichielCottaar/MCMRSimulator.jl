"""Grid acceleration over the leaf surfaces of a finite geometry tree."""
module Grids

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, InternalBoundingBox, leaf_bounding_boxes
import ..PhysicalGeometries: intersection_type, inside_indices_eltype, find_intersection
import ..PhysicalGeometries: inside_indices, isinside_single, rewrap, has_equivalent_wrapper
import ..PhysicalGeometries: to_inside_index, finalise!
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
    fully_inside_indices::Dict{SVector{N, Int}, Vector{V}}
    user_resolution::Union{Nothing, Float64}
end

function Grid(geometry::P; resolution=nothing) where {N, P<:PhysicalGeometry{N}}
    boxes = InternalBoundingBox{N}[]
    indices = Vector{intersection_type(P)}()
    _flatten_leaf_bounding_boxes!(boxes, indices, leaf_bounding_boxes(geometry), ())
    grid = IntersectionGrid(boxes; resolution, delayed=true)
    inside_type = inside_indices_eltype(P)
    fully_inside_indices = Dict{SVector{N, Int}, Vector{inside_type}}()
    Grid{N, P, eltype(indices), inside_indices_eltype(P)}(
        geometry, grid, indices, fully_inside_indices, resolution,
    )
end

function finalise!(wrapper::Grid{N}) where {N}
    finalise!(wrapper.geometry)
    boxes = InternalBoundingBox{N}[]
    _flatten_leaf_bounding_boxes!(boxes, leaf_bounding_boxes(wrapper.geometry))
    finalise!(wrapper.grid, boxes)
    empty!(wrapper.fully_inside_indices)
    merge!(
        wrapper.fully_inside_indices,
        _fully_inside_indices(wrapper.geometry, wrapper.grid, wrapper.indices),
    )
    nothing
end

function _fully_inside_indices(geometry, grid::IntersectionGrid{N}, collision_indices) where {N}
    inside_type = inside_indices_eltype(typeof(geometry))
    result = Dict{SVector{N, Int}, Vector{inside_type}}()

    inside_boxes = InternalBoundingBox{N}[]
    inside_leaf_indices = Vector{inside_type}()
    _flatten_inside_bounding_boxes!(
        inside_boxes, inside_leaf_indices, leaf_bounding_boxes(geometry; inside=true), (),
    )
    isempty(inside_boxes) && return result

    inside_grid_indices = InternalBoundingBoxes.grid_indices(
        grid.grid_bounding_box, grid.dimensions, inside_boxes,
    )
    collision_inside_indices = [
        to_inside_index(geometry, collision_index) for collision_index in collision_indices
    ]
    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    cell_size = 1 ./ grid.inv_resolution
    for voxel in CartesianIndices(Tuple(grid.dimensions))
        isempty(inside_grid_indices[voxel]) && continue
        centre = lower_bound .+ (SVector{N, Float64}(Tuple(voxel)) .- 0.5) .* cell_size
        for inside_leaf_index in inside_grid_indices[voxel]
            inside_index = inside_leaf_indices[inside_leaf_index]
            any(
                collision_index -> collision_inside_indices[collision_index] == inside_index,
                grid.indices[voxel],
            ) && continue
            if isinside_single(geometry, centre; candidate=inside_index)
                coordinate = SVector{N, Int}(Tuple(voxel))
                push!(get!(result, coordinate, inside_type[]), inside_index)
            end
        end
    end
    result
end

function _flatten_inside_bounding_boxes!(boxes, indices, box::InternalBoundingBox, path)
    push!(boxes, box)
    push!(indices, path)
    nothing
end

function _flatten_inside_bounding_boxes!(boxes, indices, children::Union{Tuple, AbstractVector}, path)
    for (index, child) in enumerate(children)
        _flatten_inside_bounding_boxes!(boxes, indices, child, (path..., index))
    end
    nothing
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

function _flatten_leaf_bounding_boxes!(boxes, box::InternalBoundingBox)
    push!(boxes, box)
    nothing
end

function _flatten_leaf_bounding_boxes!(boxes, children::Union{Tuple, AbstractVector})
    for child in children
        _flatten_leaf_bounding_boxes!(boxes, child)
    end
    nothing
end

InternalBoundingBox(geometry::Grid; kwargs...) = geometry.grid.bounding_box

function _grid_coordinate(grid::IntersectionGrid{N}, position::SVector{N, Float64}) where {N}
    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    upper_bound = InternalBoundingBoxes.upper(grid.grid_bounding_box)
    (any(position .< lower_bound) || any(position .> upper_bound)) && return nothing
    coordinate = Int.(floor.((position - lower_bound) .* grid.inv_resolution)) .+ 1
    coordinate = clamp.(coordinate, 1, SVector{N, Int}(size(grid.indices)))
    SVector{N, Int}(coordinate)
end

function inside_indices(
    wrapper::Grid{N},
    position::SVector{N, Float64},
    intersection=nothing;
    no_deproject=false,
    kwargs...,
) where {N}
    inside_type = inside_indices_eltype(typeof(wrapper.geometry))
    indices = Vector{inside_type}()
    _require_finalised(wrapper.grid)
    coordinate = _grid_coordinate(wrapper.grid, position)
    isnothing(coordinate) && return indices

    append!(indices, get(wrapper.fully_inside_indices, coordinate, inside_type[]))
    checked = Set{inside_type}()
    for collision_index in wrapper.grid.indices[coordinate...]
        inside_index = to_inside_index(wrapper.geometry, wrapper.indices[collision_index])
        isnothing(inside_index) && continue
        inside_index in checked && continue
        push!(checked, inside_index)
        isinside_single(
            wrapper.geometry, position, intersection;
            candidate=inside_index, no_deproject, kwargs...,
        ) && push!(indices, inside_index)
    end
    sort!(indices; lt=isless)
end

function _require_finalised(grid::IntersectionGrid)
    isassigned(grid.indices, lastindex(grid.indices)) || throw(ArgumentError(
        "Grid is not finalised; call finalise!(grid) before querying it",
    ))
    nothing
end

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
    _require_finalised(geometry.grid)
    iszero(destination - start) && return nothing
    start_coordinate = _grid_coordinate(geometry.grid, start)
    destination_coordinate = _grid_coordinate(geometry.grid, destination)
    if !isnothing(start_coordinate) && start_coordinate == destination_coordinate
        return _find_intersection_candidates(
            geometry, start, destination, previous_hit, inside,
            ((index, 0.0) for index in geometry.grid.indices[start_coordinate...]),
        )
    end
    _find_intersection_candidates(
        geometry, start, destination, previous_hit, inside,
        GridIterator(geometry.grid, start, destination),
    )
end

function _find_intersection_candidates(geometry, start, destination, previous_hit, inside, candidates)
    best = nothing
    for (index, entry_time) in candidates
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
