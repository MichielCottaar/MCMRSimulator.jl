module PhysicalGeometries
import StaticArrays: SVector
import ..Indices: ObstructionIndex
import ..InternalBoundingBoxes: InternalBoundingBox
import ...BoundingBoxes: BoundingBox

"""
    PhysicalGeometry{N}

Represents a physical geometry in N-dimensional space. This is an abstract type that can be extended to define specific geometries and their properties.
"""
abstract type PhysicalGeometry{N} end

function size_scale end
function geometry_mesh end

"""
    Intersection{N}(distance, normal, inside, obstruction_index, hit_gap)

Intersection between a path and a physical geometry. `distance` is normalized
to the path from start to destination and `obstruction_index` identifies the
obstruction that was hit.
"""
struct Intersection{N}
    distance::Float64
    normal::SVector{N, Float64}
    inside::Bool
    obstruction_index::ObstructionIndex
    hit_gap::Bool
end

Intersection{N}() where {N} = Intersection(Inf, zero(SVector{N, Float64}), false, ObstructionIndex(), false)
function Base.isempty(intersection::Intersection)
    intersection.distance < 0 || intersection.distance > 1
end

"""Find the first intersection of a path with `geometry`."""
function detect_intersection end

"""Return whether a geometry has an inside region."""
function has_inside end

"""Return the obstruction indices containing a position."""
function inside_indices end

"""Sample surface positions and inward normals within a bounding box."""
function surface_sampling end

"""Sample surface positions together with initialized collision states."""
function random_surface_positions end

"""Return the axis-aligned bounding box containing a physical geometry."""
function InternalBoundingBox end
function _geometry_mesh end
_mesh_result(vertices, triangles) = (vertices=vertices, triangles=triangles)
_translate_native(value, displacement) = value isa Real ? value + displacement[1] :
    value isa SVector ? value + displacement :
    value isa NamedTuple ? _mesh_result(value.vertices .+ Ref(displacement), value.triangles) :
    [_translate_native(item, displacement) for item in value]

function geometry_mesh(geometry; height=nothing, nsamples=100, bounding_box=nothing)
    bounding_box = bounding_box isa BoundingBox ? InternalBoundingBox(bounding_box) : bounding_box
    _geometry_mesh(geometry; height, nsamples, bounding_box)
end

include("grid_dispatch.jl")
include("groups.jl")
include("transformations.jl")
include("repeats.jl")
include("transparent.jl")
include("base_obstructions/base_obstructions.jl")
include("meshes.jl")
import .Transparents: Transparent, transparent_geometry

end
