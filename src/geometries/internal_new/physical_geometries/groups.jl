"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import ..PhysicalGeometries: PhysicalGeometry

abstract type GroupGeometry{N} <: PhysicalGeometry{N} end

struct GeometryVector{N, P<:PhysicalGeometry{N}} <: GroupGeometry{N}
    geometries::Vector{P}
end

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} <: GroupGeometry{N}
    geometries::P
end

GeometryVector{N}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P}(geometries)

GeometryTuple{N}(geometries::Tuple{Vararg{PhysicalGeometry{N}}}) where {N} =
    GeometryTuple{N, typeof(geometries)}(geometries)

Base.size(geometry::GeometryVector) = size(geometry.geometries)
Base.axes(geometry::GeometryVector) = axes(geometry.geometries)
Base.length(geometry::GeometryVector) = length(geometry.geometries)
Base.firstindex(geometry::GeometryVector) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryVector) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryVector, index...) = getindex(geometry.geometries, index...)
Base.setindex!(geometry::GeometryVector, value, index...) = setindex!(geometry.geometries, value, index...)
Base.iterate(geometry::GeometryVector, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{GeometryVector{N, P}}) where {N, P} = P

Base.length(geometry::GeometryTuple) = length(geometry.geometries)
Base.firstindex(geometry::GeometryTuple) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryTuple) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryTuple, index...) = getindex(geometry.geometries, index...)
Base.iterate(geometry::GeometryTuple, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{GeometryTuple{N, P}}) where {N, P} = eltype(P)
Base.Tuple(geometry::GeometryTuple) = geometry.geometries

end
