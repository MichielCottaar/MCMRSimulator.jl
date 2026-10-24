module PhysicalGeometries
import StaticArrays: SVector
import Random: rand, default_rng
import Distributions: Poisson
import ..InternalBoundingBoxes
import ..InternalBoundingBoxes: InternalBoundingBox
import ...BoundingBoxes: BoundingBox
import ..InsideViews: InsideView, child_view

"""
    PhysicalGeometry{N}

Represents a physical geometry in N-dimensional space. This is an abstract type that can be extended to define specific geometries and their properties.
"""
abstract type PhysicalGeometry{N} end

"""
    has_equivalent_wrapper(geometry1::PhysicalGeometry, geometry2::PhysicalGeometry)
    has_equivalent_wrapper(group)

Return whether two geometries have equivalent top-level wrappers, independently
of their wrapped geometries. The default is `false`, including for identical
objects whose wrapper type has no specialized comparison. Supported wrapper
comparisons should be symmetric and transitive.
`Repeat` wrappers are equivalent when their repeat-distance vectors match
exactly; their children may differ in type and content.

For a vector-like or tuple group, compare every child with the first child.
Empty groups return `false`; a singleton group compares its child with itself.
"""
function has_equivalent_wrapper end

has_equivalent_wrapper(::PhysicalGeometry, ::PhysicalGeometry) = false

"""
    rewrap(wrapper, geometry::PhysicalGeometry)

Construct the same kind of wrapper around `geometry`, preserving the wrapper's
parameters. Wrapper types supporting common-wrapper extraction implement this
alongside `has_equivalent_wrapper` and expose their child as `geometry`.
"""
function rewrap end

"""
    extract_common_wrapper(group)

Lift consecutive equivalent wrapper layers from every child of a plain
`GeometryVector` or `GeometryTuple` into wrappers around the group. After removing
each common layer, recursively extract wrappers from the resulting inner group
before reapplying the outer wrapper. Stop when `has_equivalent_wrapper` is false.
Preserve child order, group family, and wrapper order without mutating the input.

Return the original group unchanged if its children have no equivalent wrapper.
Recursion only follows inner groups created by removing common wrappers; it does
not otherwise simplify descendants.

For example, a group of `Repeat` objects with identical repeat distances becomes
a `Repeat` of the group of their children. Geometric paths change from
`(child_index, repeat_shift, ...)` to `(repeat_shift, child_index, ...)`, so cached
inside indices and previous-hit paths must be regenerated. Property paths remain
unchanged after stripping repeat indices.
"""
function extract_common_wrapper end

"""
    simplify(geometry::PhysicalGeometry)

Return an equivalent internal geometry representation with common wrappers
extracted. Simplify the children of supported groups and wrappers first, then
apply `extract_common_wrapper` to plain `GeometryVector` and `GeometryTuple`
groups. Base obstructions are returned unchanged. Every other geometry type must
provide simplification support; missing implementations raise `MethodError`.
The input is never mutated, and unchanged subtrees may be returned by identity.

A vector group accepts simplified children only if all have the same fully
concrete type. Otherwise, discard all proposed child changes and return the
original vector group; never convert it to a tuple group. Tuple groups retain
their tuple representation and may accept heterogeneous child changes.

Wrapper extraction can reorder geometric paths. Generate cached inside states
and previous-hit paths for the returned geometry rather than reusing those from
the original representation. Property paths retain their original ordering.
"""
function simplify end

struct IntersectionParams{N}
    inside::Bool
    normal::SVector{N, Float64}
    hit_gap::Bool
end

function size_scale end
function geometry_mesh end
function distance_to_surface end


"""
    find_intersection(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, previous_hit=nothing, inside=nothing, candidate=nothing) -> Optional{indices}

Finds the first intersection between `start` and `dest`.

This returns a tuple with just enough information to identify this intersection, namely:
1. an increasingly lengthy list of the `index` (int) of the closest intersection in a group or shift (`SVector{3, Int}`)
2. whether the path enters from the inside as the penultimate element
3. the distance of the intersection from `start` (0) to `dest` (1) as the final element
Each layer prepends its identifier information to the front (which can be something else).

The final two elements therefore always have the form `(..., inside, distance)`.

Return `nothing` if there is no intersection between `start` and `dest`.

`candidate` is `nothing` for unrestricted traversal, or a tuple identifying one
leaf surface in the hierarchy returned by `leaf_bounding_boxes`. Groups dispatch
directly to the selected child and consume one path component. Transformations
and transparent wrappers forward the path unchanged. Ordinary base obstructions
expect `()`, meshes expect `(triangle_index,)`, and finite cylinders expect
`(surface_index,)` (1: side, 2: first cap, 3: second cap).

Candidate selection preserves the original segment, intersection paths,
previous-hit handling, and cached `inside` state. Repeating and liminal geometries
require `candidate === nothing`; selection belongs inside their finite children.
"""
function find_intersection end
find_intersection_requires_inside(::Type{<:PhysicalGeometry}) = Val(false)

function is_liminal end
is_liminal(::Type{<:PhysicalGeometry}) = false

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
     get_intersection_params(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, indices[, isinside]) -> IntersectionParams{N}

Gets the `inside`, `normal`, and `hit_gap` properties of the intersection after it has been identified by `find_intersection`. When provided, `isinside` contains the cached obstruction indices for the current particle and is converted to child-local indices during recursion.

This should be overwritten for base obstructions, but should work as is for any composite geometries, which override `get_child` instead.
"""
function get_intersection_params(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, indices::Tuple, isinside=nothing) where {N}
    (child, remaining_indices) = get_child(geometry, indices)
    start_child = to_child_coordinates(geometry, start)
    dest_child = to_child_coordinates(geometry, dest)
    n_consumed = length(indices) - length(remaining_indices)
    prefix = n_consumed == 0 ? nothing : indices[1]
    result = get_intersection_params(
        child,
        start_child,
        dest_child,
        remaining_indices,
        child_view(get_intersection_params_requires_inside(typeof(child)), isinside, prefix),
    )
    return from_child_coordinates(geometry, result) :: IntersectionParams{N}
end

get_intersection_params_requires_inside(::Type{<:PhysicalGeometry}) = Val(false)

"""
    to_child_coordinates(geometry, position/bounding_box)

Converts the `position`/`bounding_box` from `geometry` coordinates to the coordinates from the child of `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
to_child_coordinates(::PhysicalGeometry{N}, object) where {N} = object

"""
    from_child_coordinates(geometry, position/bounding_box)

Converts the `position`/`bounding_box` from child of `geometry` coordinates to the coordinates from the `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
from_child_coordinates(::PhysicalGeometry{N}, object) where {N} = object

to_child_coordinates_normal(::PhysicalGeometry{N}, object) where {N} = object

"""Convert intersection parameters from child to parent coordinates."""
from_child_coordinates(::PhysicalGeometry, params::IntersectionParams) = params

"""Return whether a geometry has an inside region."""
function has_inside end

"""Returns whether the geometry has a single inside region."""
function has_single_inside end

"""
    inside_indices_eltype(::Type{<:PhysicalGeometry})

Return the element type of the array returned by `inside_indices`.
"""
function inside_indices_eltype end
inside_indices_eltype(::Type{Union{}}) = Union{}

"""Return the index tuple type produced by `find_intersection`."""
function intersection_type end
intersection_type(::Type{Union{}}) = Union{}

"""Return intersection index types for surfaces with nonzero density."""
function bound_intersection_type end

function _merge_types(types)
    unique_types = Set{Type}()
    for type in types
        type === Union{} && continue
        type isa Union ?
            union!(unique_types, Base.uniontypes(type)) :
            push!(unique_types, type)
    end
    Union{unique_types...}
end

"""Return the obstruction indices containing a position."""
function inside_indices end

"""Sample volume positions together with their cached inside indices."""
function volume_sampling end
function inside_sampling end
function projected_surface_area end
function normal_second_moment end
function tortuosity_tensor end
function inverse_mean_free_path end

function _volume_inside_indices(geometry::PhysicalGeometry, position; no_deproject=false)
    has_inside(typeof(geometry)) || return Tuple{}[]
    has_single_inside(typeof(geometry)) &&
        return isinside_single(geometry, position; no_deproject) ? [()] : Tuple{}[]
    inside_indices(geometry, position; no_deproject)
end

function volume_sampling(
    geometry::PhysicalGeometry{N},
    bounding_box::InternalBoundingBox{N},
    volume_density::Number,
    ; no_deproject=false,
) where {N}
    volume_density >= 0 || throw(ArgumentError("volume density must be non-negative"))
    lower_bound = InternalBoundingBoxes.lower(bounding_box)
    size = 2 .* InternalBoundingBoxes.half_size(bounding_box)
    nsamples = rand(Poisson(volume_density * prod(size)))
    positions = [SVector{N, Float64}(rand(N) .* size .+ lower_bound) for _ in 1:nsamples]
    positions, [
        _volume_inside_indices(geometry, position; no_deproject)
        for position in positions
    ]
end

function _rejection_inside_sampling(
    geometry::PhysicalGeometry{N},
    density::Real,
    rng;
    no_deproject=false,
) where {N}
    bounding_box = InternalBoundingBox(geometry; no_deproject)
    lower_bound = InternalBoundingBoxes.lower(bounding_box)
    size = 2 .* InternalBoundingBoxes.half_size(bounding_box)
    box_volume = prod(size)
    nsamples = rand(rng, Poisson(density * box_volume))
    [
        position for position in (
            SVector{N, Float64}(rand(rng, N) .* size .+ lower_bound)
            for _ in 1:nsamples
        )
        if isinside_single(geometry, position; no_deproject)
    ]
end

function _inside_sampling_proposals(
    geometry::PhysicalGeometry,
    density::Real,
    rng;
    no_deproject=false,
)
    has_single_inside(typeof(geometry)) || throw(ArgumentError(
        "inside_sampling is not implemented for $(typeof(geometry))",
    ))
    _rejection_inside_sampling(geometry, density, rng; no_deproject)
end

"""Sample points uniformly within a geometry, returning positions and inside indices."""
function inside_sampling(
    geometry::PhysicalGeometry{N},
    density::Real;
    rng=default_rng(),
    no_deproject=false,
) where {N}
    isfinite(density) && density >= 0 ||
        throw(ArgumentError("inside-sampling density must be finite and nonnegative"))
    proposals = _inside_sampling_proposals(
        geometry, density, rng; no_deproject,
    )
    inside_type = inside_indices_eltype(typeof(geometry))
    positions = SVector{N, Float64}[]
    indices = Vector{inside_type}[]
    for position in proposals
        inside = has_single_inside(typeof(geometry)) ?
            (isinside_single(geometry, position; no_deproject) ? [()] : Tuple[]) :
            inside_indices(geometry, position; no_deproject)
        isempty(inside) && continue
        rand(rng) < inv(length(inside)) || continue
        push!(positions, position)
        push!(indices, inside)
    end
    positions, indices
end

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
    isempty(indices) && return ()
    child, child_indices = get_child(geometry, indices)
    cleaned = to_property_index(child, child_indices)
    nremoved = length(indices) - length(child_indices)
    return (indices[1:nremoved]..., cleaned...)
end

"""Convert a surface intersection path into an inside-region path."""
function to_inside_index(geometry::PhysicalGeometry, indices)
    has_inside(typeof(geometry)) || return nothing
    has_single_inside(typeof(geometry)) && return ()
    child, child_indices = get_child(geometry, indices)
    converted = to_inside_index(child, child_indices)
    isnothing(converted) && return nothing
    n_consumed = length(indices) - length(child_indices)
    (indices[1:n_consumed]..., converted...)
end

"""Sample surface positions together with initialized collision states."""
function random_surface_positions end

"""Return the axis-aligned bounding box containing a physical geometry."""
function InternalBoundingBox end

"""
    leaf_bounding_boxes(geometry::PhysicalGeometry)

Return conservative bounding boxes for the leaf surfaces of `geometry`, all in
the geometry's coordinate frame. The result mirrors the intersection-index
hierarchy: vector and tuple groups return vectors and tuples of child results,
while transformations and transparent wrappers preserve their child's structure.
Empty groups return empty collections.

A base obstruction returns its `InternalBoundingBox` directly. `FiniteCylinder`
returns a tuple for its side, first cap, and second cap, in surface-index order.
Meshes return a vector of triangle boxes, including gap-closing triangles.

Periodic geometries and dimension-reducing transformations have no finite leaf
bounding boxes in their outer frame. Query their finite child geometry instead.
Likewise, query individual cell geometries for liminal geometry.
"""
function leaf_bounding_boxes end

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

function Base.show(io::IO, geometry::PhysicalGeometry)
    show(io, typeof(geometry))
    show(io, "(...)")
end

include("grid_dispatch.jl")
include("groups.jl")
include("transformations.jl")
include("repeats.jl")
include("transparent.jl")
include("base_obstructions/base_obstructions.jl")
include("meshes.jl")
include("grids.jl")
include("monte_carlo_measures.jl")
import .MonteCarloMeasures: SurfaceEstimate, contains_repeat, estimate_volume, estimate_surface
include("liminal_geometries.jl")
import .Transparents: Transparent, IgnoreOverlapping, transparent_geometry
import .LiminalGeometries: FixedLiminalGeometry

end
