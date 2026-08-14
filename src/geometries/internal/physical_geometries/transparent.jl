module Transparents

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox

abstract type Transparent{N, P <: PhysicalGeometry{N}} <: PhysicalGeometry{N} end

transparent_geometry(wrapper::Transparent) = getfield(wrapper, :geometry)

has_inside(::Type{<:Transparent{N, P}}) where {N, P} = has_inside(P)

InternalBoundingBox(wrapper::Transparent) = InternalBoundingBox(transparent_geometry(wrapper))

inside_indices(
    wrapper::Transparent{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N} = inside_indices(transparent_geometry(wrapper), position, intersection)

function detect_intersection(
    wrapper::Transparent{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{N}=Intersection{N}(),
) where {N}
    detect_intersection(transparent_geometry(wrapper), start, destination, previous_hit)
end

end
