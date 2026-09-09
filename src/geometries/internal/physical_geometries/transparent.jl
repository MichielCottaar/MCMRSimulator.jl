module Transparents

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, IntersectionParams, child_type, find_intersection, get_child, get_intersection_params, get_intersection_params_requires_inside, has_inside, has_single_inside, inside_indices_eltype, intersection_type, bound_intersection_type, isinside_single, inside_indices, InternalBoundingBox, size_scale
import ..PhysicalGeometries: random_surface_positions, _geometry_mesh
import ...Properties: GeometryProperties
import ...InsideViews: has_other

abstract type Transparent{N, P <: PhysicalGeometry{N}} <: PhysicalGeometry{N} end

get_intersection_params_requires_inside(::Type{<:Transparent{N, P}}) where {N, P} =
    get_intersection_params_requires_inside(P)

"""Treat intersections inside another obstruction as gaps."""
struct IgnoreOverlapping{N, P <: PhysicalGeometry{N}} <: Transparent{N, P}
    geometry::P
end

get_intersection_params_requires_inside(::Type{<:IgnoreOverlapping}) = Val(true)

function Base.show(io::IO, ::Type{T}) where {N, P, T <: Transparent{N, P}}
    print(io, nameof(T), "{")
    show(io, P)
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

has_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_inside(P)
has_single_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_single_inside(P)
InternalBoundingBox(wrapper::Transparent) = InternalBoundingBox(transparent_geometry(wrapper))

inside_indices(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection=nothing,
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

function get_intersection_params(
    wrapper::IgnoreOverlapping{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    indices::Tuple,
    isinside=nothing,
) where {N}
    result = get_intersection_params(
        transparent_geometry(wrapper),
        start,
        destination,
        indices,
        isinside,
    )
    obstruction_indices = indices[1:(end - 2)]
    overlapping = has_other(isinside, obstruction_indices)
    IntersectionParams{N}(result.inside, result.normal, result.hit_gap || overlapping)
end

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
