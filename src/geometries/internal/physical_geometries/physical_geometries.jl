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
    find_intersection(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, previous_hit=nothing) -> Optional{indices}

Finds the first intersection between `start` and `dest`.

This returns a tuple with just enough information to identify this intersection, namely:
1. an increasingly lengthy list of the `index` (int) of the closest intersection in a group or shift (`SVector{3, Int}`)
2. whether the path enters from the inside as the penultimate element
3. the distance of the intersection from `start` (0) to `dest` (1) as the final element
Each layer prepends its identifier information to the front (which can be something else).

The final two elements therefore always have the form `(..., inside, distance)`.

Return `nothing` if there is no intersection between `start` and `dest`.
"""
function find_intersection end

"""
    get_child(geometry::PhysicalGeometry, indices) -> (PhysicalGeometry, remaining_indices)

Gets the child of `geometry` corresponding to the `indices` identified by `find_intersection`.

This is a helper function for `get_intersection_params`.
It should be overwritten by any wrapper obstructions (transformations, groups, repeats, etc.).
`remaining_indices` should be a tuple with any indices remaining after `get_child`.

It will not be called for base obstructions, which override `get_intersection_params` instead.
"""
function get_child end

"""
    child_type(::Type{<:PhysicalGeometry})

Return the type of the child geometry represented by a physical geometry type.
"""
function child_type end


"""
    get_intersection_params(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, indices) -> (inside=true/false, normal=SVector{N, Float64}, hit_gap=true/false)

Gets the `inside`, `normal`, and `hit_gap` properties of the intersection after it has been identified by `find_intersection`.

This should be overwritten for base obstructions, but should work as is for any composite geometries, which override `get_child` instead.
"""
function get_intersection_params(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, indices::Tuple) where {N}
    (child, remaining_indices) = get_child(geometry, indices)
    start_child = to_child_coordinates(geometry, start)
    dest_child = to_child_coordinates(geometry, dest)
    result = get_intersection_params(child, start_child, dest_child, remaining_indices)
    return (
        inside=result.inside,
        normal=from_child_coordinates_normal(geometry, result.normal),
        hit_gap=result.hit_gap,
    )
end

"""
    to_child_coordinates(geometry, position/bounding_box)

Converts the `position`/`bounding_box` from `geometry` coordinates to the coordinates from the child of `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
to_child_coordinates(::PhysicalGeometry{N}, object) = object

"""
    from_child_coordinates(geometry, position/bounding_box)

Converts the `position`/`bounding_box` from child of `geometry` coordinates to the coordinates from the `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
from_child_coordinates(::PhysicalGeometry{N}, object) = object

"""
    to_child_coordinates_normal(geometry, normal)

Converts the `normal` from `geometry` coordinates to the coordinates from the child of `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
to_child_coordinates_normal(::PhysicalGeometry{N}, object) = object

"""
    from_child_coordinates_normal(geometry, normal)

Converts the `normal` from child of `geometry` coordinates to the coordinates from the `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
from_child_coordinates_normal(::PhysicalGeometry{N}, object) = object

"""Return whether a geometry has an inside region."""
function has_inside end

"""Returns whether the geometry has a single inside region."""
function has_single_inside end

"""
    inside_indices_eltype(::Type{<:PhysicalGeometry})

Return the element type of the array returned by `inside_indices`.
"""
function inside_indices_eltype end

"""Return the obstruction indices containing a position."""
function inside_indices end

"""Whether the geometry is within the single inside."""
function isinside_single end

"""
    to_property_index(geometry, indices) -> indices

Returns a version of `indices` that is appropriate for property sampling (i.e., does not include the `repeat` index).
"""
function to_property_index(geometry::PhysicalGeometry, indices)
    if isnothing(indices)
        return nothing
    end
    child, child_indices = get_child(geometry, indices)
    cleaned = to_property_index(child, child_indices)
    nremoved = length(indices) - length(child_indices)
    return (indices[1:nremoved]..., cleaned...)
end

"""Return the obstruction-index depth of intersections from a geometry."""
function intersection_index_length end

"""Return the obstruction-index depth of inside results from a geometry."""
function inside_index_length end

"""Whether all inside-index results from a geometry have one depth."""
function all_equal_inside_depth end

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

include("intersections.jl")
import .Intersections: Intersection, flip
include("grid_dispatch.jl")
include("groups.jl")
include("transformations.jl")
include("repeats.jl")
include("transparent.jl")
include("base_obstructions/base_obstructions.jl")
include("meshes.jl")
import .Transparents: Transparent, transparent_geometry

end
