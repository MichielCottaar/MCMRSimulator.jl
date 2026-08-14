module SizeScales

"""Dynamic size-scale calculation for fixed physical geometries."""

import ..Internal: size_scale, FixedGeometry
import ..PhysicalGeometries
import ..PhysicalGeometries: PhysicalGeometry, Transparent, transparent_geometry
import ..PhysicalGeometries.BaseObstructions: InfiniteWall, Round
import ..PhysicalGeometries.Groups: GeometryVectorLike, GeometryTuple
import ..PhysicalGeometries.Transformations: Shift, Scale, Rotate
import ..PhysicalGeometries.Repeats: Repeat
import ..PhysicalGeometries.Meshes: Mesh

struct SizeScaleOverride{N, P <: PhysicalGeometry{N}} <: Transparent{N, P}
    geometry::P
    size_scale::Float64
end

function SizeScaleOverride(geometry::P, size_scale::Real) where {N, P <: PhysicalGeometry{N}}
    size_scale > 0 || throw(ArgumentError("size_scale must be positive"))
    SizeScaleOverride{N, P}(geometry, Float64(size_scale))
end

size_scale(::InfiniteWall) = Inf
size_scale(round::Round) = round.radius

function size_scale(mesh::Mesh)
    ntriangles = mesh.first_index_of_gap - 1
    ntriangles <= 1 && return Inf
    1 / (2 * PhysicalGeometries.Meshes.curvature(mesh))
end

function size_scale(geometry::GeometryVectorLike)
    isempty(geometry) && return Inf
    minimum(size_scale, geometry)
end

function size_scale(geometry::GeometryTuple)
    isempty(geometry) && return Inf
    minimum(size_scale, geometry)
end

size_scale(geometry::Shift) = size_scale(geometry.geometry)
size_scale(geometry::Rotate) = size_scale(geometry.geometry)
size_scale(geometry::Scale) = geometry.scale * size_scale(geometry.geometry)

function size_scale(geometry::Repeat)
    min(size_scale(geometry.geometry), minimum(geometry.repeats))
end

size_scale(wrapper::Transparent) = size_scale(transparent_geometry(wrapper))
size_scale(wrapper::SizeScaleOverride) = wrapper.size_scale

size_scale(geometry::FixedGeometry) = size_scale(geometry.geometry)

end
