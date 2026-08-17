module Transparents

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, has_single_inside, isinside_single, inside_indices, InternalBoundingBox, size_scale
import ..PhysicalGeometries: intersection_index_length, inside_index_length, contains_geometry_tuple
import ..PhysicalGeometries: random_surface_positions, _geometry_mesh
import ...Properties: GeometryProperties

abstract type Transparent{N, P <: PhysicalGeometry{N}} <: PhysicalGeometry{N} end

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
contains_geometry_tuple(::Type{<:Transparent{N, P}}) where {N, P} = contains_geometry_tuple(P)

InternalBoundingBox(wrapper::Transparent) = InternalBoundingBox(transparent_geometry(wrapper))

inside_indices(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N} = inside_indices(transparent_geometry(wrapper), position, intersection)

isinside_single(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N} = isinside_single(transparent_geometry(wrapper), position, intersection)

function detect_intersection(
    wrapper::Transparent{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{N}=Intersection{N}(),
) where {N}
    detect_intersection(transparent_geometry(wrapper), start, destination, previous_hit)
end

size_scale(wrapper::SizeScaleOverride) = wrapper.size_scale
size_scale(wrapper::Transparent) = size_scale(transparent_geometry(wrapper))

random_surface_positions(wrapper::Transparent{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where N =
    random_surface_positions(transparent_geometry(wrapper), density, bounding_box, scale_density)
_geometry_mesh(wrapper::Transparent; kwargs...) = _geometry_mesh(transparent_geometry(wrapper); kwargs...)

end
