"""Grid acceleration over the leaf surfaces of a finite geometry tree."""
module Grids

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, InternalBoundingBox, leaf_bounding_boxes
import ..PhysicalGeometries: find_intersection, inside_indices, inside_indices_eltype, finalise!
import ..PhysicalGeometries: filter_bounding_box, rewrap, has_equivalent_wrapper
import ..PhysicalGeometries.Groups: inside_indices_for_any_type
import ..Transparents: Transparent
import ..GridDispatch: IntersectionGrid, GridIterator
import ...InternalBoundingBoxes
import .....DebugCounters

function __init__()
    DebugCounters.register!(Grid, [
        (:ncells, "Number of grid cells between start and end of trajectory in single steps; 0 means that the start and end were identical, so no collision checks were required. This is only counted if start and destination are both inside of the grid", DebugCounters.FrequencyCounter(Int)),
        (:cell_occupancy, "Number of leaf surfaces in each cell geometry", DebugCounters.FrequencyCounter(Int)),
        (:cells_visited, "Number of grid cells actually visited before a collision was found. Empty cells will not be counted.", DebugCounters.FrequencyTable(Int)),
        (:ncandidates, "Number of non-empty cell geometries checked by the collision detector", DebugCounters.FrequencyTable(Int)),
        (:duplicate_candidates, "Legacy candidate counter; geometry payload traversal records zero duplicates", DebugCounters.FrequencyTable(Int)),
        (:start_out_of_grid, "Is the start position outside of the grid", DebugCounters.FrequencyTable(Bool)),
        (:dest_out_of_grid, "Is the destination position outside of the grid", DebugCounters.FrequencyTable(Bool)),
    ])
end

"""
    Grid(geometry::PhysicalGeometry; resolution=nothing, debug=false)

Wrap a finite geometry in an intersection grid over its leaf bounding boxes.
`resolution` specifies the cell size; `nothing` uses the `IntersectionGrid`
estimate. The requested value is stored in `user_resolution` and reused when
rewrapping a simplified child. Leaf indices follow the hierarchy returned by
`leaf_bounding_boxes`, so intersection paths and properties retain their original
geometry indices. Each cell stores a filtered geometry subset, or `nothing` if
the cell does not overlap the geometry.

The wrapper expects an unrestricted incoming candidate and supplies each visited
cell's filtered geometry to its normal query methods. For periodic geometries,
place the grid inside `Repeat`. Set `debug=true` to collect per-instance debug
metrics.
"""
struct Grid{N, P<:PhysicalGeometry{N}, D} <: Transparent{N, P}
    geometry::P
    grid::IntersectionGrid{N}
    user_resolution::Union{Nothing, Float64}
end

function Base.show(io::IO, ::Type{Grid{N, P, D}}) where {N, P, D}
    print(io, "Grid{")
    show(io, P)
    print(io, "}")
end

function Grid(geometry::P; resolution=nothing, debug=false) where {N, P<:PhysicalGeometry{N}}
    debug isa Bool || throw(ArgumentError("debug must be a Bool"))
    boxes = InternalBoundingBox{N}[]
    _flatten_leaf_bounding_boxes!(boxes, leaf_bounding_boxes(geometry))
    grid = IntersectionGrid(boxes; resolution)
    Grid{N, P, debug}(geometry, grid, resolution)
end

function finalise!(wrapper::Grid{N}) where {N}
    finalise!(wrapper.geometry)
    finalise!(wrapper.grid, wrapper.geometry)
    if DebugCounters.debug_mode(wrapper)
        for geometry in wrapper.grid.geometries
            DebugCounters.record!(wrapper, :cell_occupancy, _cell_leaf_count(geometry))
        end
    end
    nothing
end

_cell_leaf_count(::Nothing) = 0

_cell_leaf_count(box::InternalBoundingBox) = 1

_cell_leaf_count(boxes::Union{Tuple, AbstractVector}) =
    sum(_cell_leaf_count, boxes; init=0)

_cell_leaf_count(geometry::PhysicalGeometry) = _cell_leaf_count(leaf_bounding_boxes(geometry))

DebugCounters.debug_mode(::Grid{N, P, D}) where {N, P, D} = D
rewrap(wrapper::Grid{N}, geometry::PhysicalGeometry{N}) where {N} =
    Grid(geometry; resolution=wrapper.user_resolution, debug=DebugCounters.debug_mode(wrapper))
has_equivalent_wrapper(first::Grid{N}, second::Grid{N}) where {N} =
    first.user_resolution == second.user_resolution

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
    coordinate = clamp.(coordinate, 1, grid.dimensions)
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
    cell_geometry = wrapper.grid.geometries[coordinate...]
    isnothing(cell_geometry) && return indices
    append!(indices, inside_indices_for_any_type(
        cell_geometry, position, intersection; no_deproject,
    ))
    sort!(indices; lt=isless)
end

function _require_finalised(grid::IntersectionGrid)
    isassigned(grid.geometries, lastindex(grid.geometries)) || throw(ArgumentError(
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
    if iszero(destination - start)
        DebugCounters.record!(geometry, :ncells, 0)
        DebugCounters.record!(geometry, :ncandidates, 0)
        DebugCounters.record!(geometry, :duplicate_candidates, 0)
        DebugCounters.record!(geometry, :cells_visited, 0)
        return nothing
    end
    start_coordinate = _grid_coordinate(geometry.grid, start)
    destination_coordinate = _grid_coordinate(geometry.grid, destination)
    if DebugCounters.debug_mode(geometry) 
        DebugCounters.record!(geometry, :start_out_of_grid, isnothing(start_coordinate))
        DebugCounters.record!(geometry, :dest_out_of_grid, isnothing(destination_coordinate))
        if !isnothing(start_coordinate) && !isnothing(destination_coordinate)
            DebugCounters.record!(geometry, :ncells, sum(abs.(destination_coordinate - start_coordinate)) + 1)
        end
    end
    if !isnothing(start_coordinate) && start_coordinate == destination_coordinate
        cell_geometry = geometry.grid.geometries[start_coordinate...]
        return isnothing(cell_geometry) ? nothing : _find_intersection_candidates(
            geometry, start, destination, previous_hit, inside, ((cell_geometry, 0.0),),
        )
    end
    _find_intersection_candidates(
        geometry, start, destination, previous_hit, inside,
        GridIterator(geometry.grid, start, destination),
    )
end

function _find_intersection_candidates(geometry, start, destination, previous_hit, inside, candidates)
    best = nothing
    if DebugCounters.debug_mode(geometry)
        ncandidates = 0
        cells_visited = 0
        prev_entry_time = -Inf
    end
    for (cell_geometry, entry_time) in candidates
        !isnothing(best) && best[end] < entry_time && break
        if DebugCounters.debug_mode(geometry)
            ncandidates += 1
            if prev_entry_time != entry_time
                cells_visited += 1
                prev_entry_time = entry_time
            end
        end
        hit = find_intersection(cell_geometry, start, destination, previous_hit, inside)
        if !isnothing(hit) && _earlier_hit(hit, best)
            best = hit
        end
    end
    if DebugCounters.debug_mode(geometry)
        DebugCounters.record!(geometry, :ncandidates, ncandidates)
        DebugCounters.record!(geometry, :duplicate_candidates, 0)
        DebugCounters.record!(geometry, :cells_visited, cells_visited)
    end
    best
end

end
