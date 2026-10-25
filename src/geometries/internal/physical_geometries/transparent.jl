module Transparents

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, IntersectionParams, child_type, find_intersection, find_intersection_requires_inside, get_child, get_intersection_params, get_intersection_params_requires_inside, has_inside, has_single_inside, inside_indices_eltype, intersection_type, bound_intersection_type, isinside_single, inside_indices, InternalBoundingBox, size_scale, to_inside_index, to_property_index
import ..PhysicalGeometries: random_surface_positions, _geometry_mesh
import ..PhysicalGeometries: _inside_sampling_proposals
import ..PhysicalGeometries: leaf_bounding_boxes
import ..PhysicalGeometries: distance_to_surface
import ..PhysicalGeometries: simplify, rewrap, rewrap_extracted_wrapper, has_equivalent_wrapper
import ...Properties: GeometryProperties
import ...InsideViews: has_other

abstract type Transparent{N, P <: PhysicalGeometry{N}} <: PhysicalGeometry{N} end

_inside_sampling_proposals(wrapper::Transparent, density::Real, rng; kwargs...) =
    _inside_sampling_proposals(
        transparent_geometry(wrapper), density, rng; kwargs...,
    )

get_intersection_params_requires_inside(::Type{<:Transparent{N, P}}) where {N, P} =
    get_intersection_params_requires_inside(P)

"""Treat intersections inside another obstruction as gaps.

`K` is the number of vector/tuple group levels this wrapper has been extracted
through. For `K > 0`, gap detection only considers other inside regions sharing
the same first `K` non-repeat indices, preserving the scope of the original
wrapper placement.
"""
struct IgnoreOverlapping{N, P <: PhysicalGeometry{N}, K} <: Transparent{N, P}
    geometry::P
end

IgnoreOverlapping(geometry::P) where {N, P<:PhysicalGeometry{N}} =
    IgnoreOverlapping{N, P, 0}(geometry)

has_equivalent_wrapper(::IgnoreOverlapping{N, P, K}, ::IgnoreOverlapping{N, Q, K}) where {N, P, Q, K} = true

get_intersection_params_requires_inside(::Type{<:IgnoreOverlapping}) = Val(true)

function Base.show(io::IO, ::Type{T}) where {T <: Transparent}
    print(io, nameof(T))
    T isa UnionAll && return
    print(io, "{")
    show(io, T.parameters[2])
    print(io, "}")
end

child_type(::Type{<:Transparent{N, P}}) where {N, P} = P
inside_indices_eltype(::Type{<:Transparent{N, P}}) where {N, P} = inside_indices_eltype(P)
intersection_type(::Type{<:Transparent{N, P}}) where {N, P} = intersection_type(P)
bound_intersection_type(wrapper::Transparent, density) =
    bound_intersection_type(transparent_geometry(wrapper), density)

struct SizeScaleOverride{N, P <: PhysicalGeometry{N}} <: Transparent{N, P}
    geometry::P
    size_scale::Float64
end

function SizeScaleOverride(geometry::P, size_scale::Real) where {N, P <: PhysicalGeometry{N}}
    size_scale > 0 || throw(ArgumentError("size_scale must be positive"))
    SizeScaleOverride{N, P}(geometry, Float64(size_scale))
end

transparent_geometry(wrapper::Transparent) = getfield(wrapper, :geometry)

function simplify(wrapper::Transparent, pass::Val)
    geometry = simplify(transparent_geometry(wrapper), pass)
    geometry === transparent_geometry(wrapper) ? wrapper : rewrap(wrapper, geometry)
end

rewrap(wrapper::IgnoreOverlapping{N, P, K}, geometry::PhysicalGeometry{N}) where {N, P, K} =
    IgnoreOverlapping{N, typeof(geometry), K}(geometry)

function rewrap_extracted_wrapper(
    reference::IgnoreOverlapping{N, P, K},
    references::Union{AbstractVector, Tuple},
    geometry::PhysicalGeometry{N},
) where {N, P, K}
    all(wrapper -> has_equivalent_wrapper(reference, wrapper), references) ||
        throw(ArgumentError("cannot extract IgnoreOverlapping wrappers with different scopes"))
    IgnoreOverlapping{N, typeof(geometry), K + 1}(geometry)
end
rewrap(wrapper::SizeScaleOverride{N}, geometry::PhysicalGeometry{N}) where {N} =
    SizeScaleOverride(geometry, wrapper.size_scale)
has_equivalent_wrapper(first::SizeScaleOverride{N}, second::SizeScaleOverride{N}) where {N} =
    first.size_scale == second.size_scale

leaf_bounding_boxes(wrapper::Transparent) = leaf_bounding_boxes(transparent_geometry(wrapper))
distance_to_surface(wrapper::Transparent, position) =
    distance_to_surface(transparent_geometry(wrapper), position)

has_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_single_inside(P)
InternalBoundingBox(wrapper::Transparent; kwargs...) =
    InternalBoundingBox(transparent_geometry(wrapper); kwargs...)

inside_indices(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection=nothing,
    ; no_deproject=false,
    kwargs...,
) where {N} = inside_indices(
    transparent_geometry(wrapper), position, intersection; no_deproject, kwargs...
)

isinside_single(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    previous_intersection=nothing,
    ; no_deproject=false,
    kwargs...,
) where {N} = isinside_single(
    transparent_geometry(wrapper), position, previous_intersection;
    no_deproject, kwargs...
)

function find_intersection(
    wrapper::Transparent{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
    inside=nothing,
    candidate::Union{Nothing, Tuple}=nothing,
) where {N}
    find_intersection(transparent_geometry(wrapper), start, destination, previous_hit, inside, candidate)
end

find_intersection_requires_inside(::Type{<:Transparent{N, P}}) where {N, P} =
    find_intersection_requires_inside(P)

get_child(wrapper::Transparent, indices) = (transparent_geometry(wrapper), indices)

function get_intersection_params(
    wrapper::IgnoreOverlapping{N, P, K},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    indices::Tuple,
    isinside=nothing,
) where {N, P, K}
    result = get_intersection_params(
        transparent_geometry(wrapper),
        start,
        destination,
        indices,
        isinside,
    )
    collision_indices = indices[1:(end - 2)]
    obstruction_indices = to_inside_index(transparent_geometry(wrapper), collision_indices)
    overlapping = !isnothing(obstruction_indices) && _has_scoped_other(
        Val(K), transparent_geometry(wrapper), isinside, obstruction_indices,
    )
    IntersectionParams{N}(result.inside, result.normal, result.hit_gap || overlapping)
end

size_scale(wrapper::SizeScaleOverride) = wrapper.size_scale
size_scale(wrapper::Transparent) = size_scale(transparent_geometry(wrapper))

function random_surface_positions(
    wrapper::Transparent{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density;
    include_gap=false,
    kwargs...,
) where {N}
    random_surface_positions(
        transparent_geometry(wrapper), density, bounding_box, scale_density;
        include_gap, kwargs...,
    )
end

function random_surface_positions(
    wrapper::IgnoreOverlapping{N, P, K},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density;
    include_gap=false,
    kwargs...,
) where {N, P, K}
    geometry = transparent_geometry(wrapper)
    positions, indices = random_surface_positions(
        geometry, density, bounding_box, scale_density; include_gap, kwargs...,
    )
    inside_indices_at(position) = has_single_inside(typeof(geometry)) ?
        (isinside_single(geometry, position) ? [()] : Tuple{}[]) :
        inside_indices(geometry, position)
    keep = [if K == 0
        # Retain the legacy surface-path comparison when no group was lifted.
        all(inside_index == index[1:(end - 1)] for inside_index in inside_indices_at(position))
    else
        obstruction_indices = to_inside_index(geometry, index[1:(end - 1)])
        isnothing(obstruction_indices) || !_has_scoped_other(
            Val(K), geometry, inside_indices_at(position), obstruction_indices,
        )
    end for (position, index) in zip(positions, indices)]
    positions[keep], indices[keep]
end

# Preserve the original, allocation-free comparison when no group wrapper has
# been extracted. Lifted wrappers compare prefixes after removing repeat indices.
_has_scoped_other(::Val{0}, geometry, inside, obstruction_indices) =
    has_other(inside, obstruction_indices)

function _has_scoped_other(::Val{K}, geometry, inside, obstruction_indices) where {K}
    obstruction_property_index = to_property_index(geometry, obstruction_indices)
    (isnothing(obstruction_property_index) || length(obstruction_property_index) < K) && return false
    obstruction_prefix = _index_prefix(obstruction_property_index, K)
    for inside_index in inside
        property_index = to_property_index(geometry, inside_index)
        (isnothing(property_index) || length(property_index) < K) && continue
        property_index == obstruction_property_index && continue
        _index_prefix(property_index, K) == obstruction_prefix && return true
    end
    false
end

_index_prefix(index::Tuple, count::Int) = index[1:min(count, length(index))]

_geometry_mesh(wrapper::Transparent; kwargs...) = _geometry_mesh(transparent_geometry(wrapper); kwargs...)

end
