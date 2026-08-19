module Transparents

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, child_type, find_intersection, get_child, has_inside, has_single_inside, inside_indices_eltype, isinside_single, inside_indices, InternalBoundingBox, size_scale
import ..PhysicalGeometries: intersection_index_length, inside_index_length, all_equal_inside_depth
import ..PhysicalGeometries: random_surface_positions, _geometry_mesh
import ...Properties: GeometryProperties

abstract type Transparent{N, P <: PhysicalGeometry{N}} <: PhysicalGeometry{N} end

child_type(::Type{<:Transparent{N, P}}) where {N, P} = P
inside_indices_eltype(::Type{<:Transparent{N, P}}) where {N, P} = inside_indices_eltype(P)

struct SizeScaleOverride{N, P <: PhysicalGeometry{N}} <: Transparent{N, P}
    geometry::P
    size_scale::Float64
end

function SizeScaleOverride(geometry::P, size_scale::Real) where {N, P <: PhysicalGeometry{N}}
    size_scale > 0 || throw(ArgumentError("size_scale must be positive"))
    SizeScaleOverride{N, P}(geometry, Float64(size_scale))
end

transparent_geometry(wrapper::Transparent) = getfield(wrapper, :geometry)

has_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_single_inside(P)
intersection_index_length(::Type{<:Transparent{N, P}}) where {N, P} = intersection_index_length(P)
inside_index_length(::Type{<:Transparent{N, P}}) where {N, P} = inside_index_length(P)
all_equal_inside_depth(::Type{<:Transparent{N, P}}) where {N, P} = all_equal_inside_depth(P)

InternalBoundingBox(wrapper::Transparent) = InternalBoundingBox(transparent_geometry(wrapper))

inside_indices(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection::Intersection{3}=Intersection{3}(),
) where {N} = inside_indices(transparent_geometry(wrapper), position, intersection)

isinside_single(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    previous_intersection=nothing,
) where {N} = isinside_single(transparent_geometry(wrapper), position, previous_intersection)

function find_intersection(
    wrapper::Transparent{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
) where {N}
    find_intersection(transparent_geometry(wrapper), start, destination, previous_hit)
end

get_child(wrapper::Transparent, indices) = (transparent_geometry(wrapper), indices)

size_scale(wrapper::SizeScaleOverride) = wrapper.size_scale
size_scale(wrapper::Transparent) = size_scale(transparent_geometry(wrapper))

function random_surface_positions(
    wrapper::Transparent{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    random_surface_positions(transparent_geometry(wrapper), density, bounding_box, scale_density)
end
_geometry_mesh(wrapper::Transparent; kwargs...) = _geometry_mesh(transparent_geometry(wrapper); kwargs...)

end
