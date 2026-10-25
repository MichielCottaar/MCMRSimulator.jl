"""Grid acceleration over the leaf surfaces of a finite geometry tree."""
module Grids

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, InternalBoundingBox, leaf_bounding_boxes
import ..PhysicalGeometries: intersection_type, find_intersection, rewrap, has_equivalent_wrapper
import ..Transparents: Transparent
import ..GridDispatch: IntersectionGrid, GridIterator

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
struct Grid{N, P<:PhysicalGeometry{N}, I} <: Transparent{N, P}
    geometry::P
    grid::IntersectionGrid{N}
    indices::Vector{I}
    user_resolution::Union{Nothing, Float64}
end

function Grid(geometry::P; resolution=nothing) where {N, P<:PhysicalGeometry{N}}
    boxes = InternalBoundingBox{N}[]
    indices = Vector{intersection_type(P)}()
    _flatten_leaf_bounding_boxes!(boxes, indices, leaf_bounding_boxes(geometry), ())
    grid = IntersectionGrid(boxes; resolution)
    Grid{N, P, eltype(indices)}(geometry, grid, indices, resolution)
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
