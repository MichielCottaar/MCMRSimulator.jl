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
    child_geometry(geometry::PhysicalGeometry)

Return the single wrapped child geometry, or throw if the geometry has no single child.
"""
function child_geometry end

function child_geometry(geometry::PhysicalGeometry)
    throw(ArgumentError("$(typeof(geometry)) does not have a single child geometry"))
end

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
    rewrap(group, children)
    rewrap(wrapper, geometry::PhysicalGeometry, pass::Val)

Construct the same kind of wrapper around `geometry`, preserving the wrapper's
parameters. Wrapper types supporting common-wrapper extraction implement this
alongside `has_equivalent_wrapper` and expose their child as `geometry`.
For vector and tuple groups, reconstruct the same group family from the supplied
child collection, rebuilding acceleration data where applicable.
"""
function rewrap end
rewrap(wrapper, geometry::PhysicalGeometry, ::Val) = rewrap(wrapper, geometry)

"""
    rewrap_extracted_wrapper(wrapper, references, geometry::PhysicalGeometry)

Rebuild an extracted wrapper around `geometry`, using source `references` when
the wrapper needs to adjust its indexing scope.
"""
rewrap_extracted_wrapper(wrapper, references, geometry::PhysicalGeometry) =
    rewrap(wrapper, geometry)

"""
    simplify(geometry::PhysicalGeometry)
    simplify(geometry::PhysicalGeometry, pass::Val)

Return an equivalent internal geometry representation after applying the
simplification pass identified by `Val(:symbol)`. Default implementations forward
the same pass to their children and reconstruct through `rewrap` when needed.
Base obstructions are returned unchanged. Pass-specific methods can add local
rewrites; an otherwise unrecognized pass just traverses supported nodes.
Geometry types without traversal support raise `MethodError`. The input is never
mutated, and unchanged subtrees may be returned by identity.

A vector group accepts simplified children only if all have the same fully
concrete type. Otherwise, discard all proposed child changes and use the original
vector group for any local rewrite; never convert it to a tuple group. Rejecting
child changes does not skip the local rewrite. Tuple groups retain their tuple
representation and may accept heterogeneous child changes.
Bounding-box vector groups follow the same child-type policy, rebuilding their
boxes when children change.

`Val(:extract_common_wrapper)` performs one bottom-up pass, extracting at most one
common wrapper layer from each visited plain vector or tuple group. For example,
a group of `Repeat` objects with identical periods becomes a `Repeat` of the group
of their children. The newly constructed inner group is left for a subsequent
pass. Child and wrapper order are preserved.

Without a pass argument, repeatedly apply `Val(:extract_common_wrapper)` until
the returned representation is unchanged according to `===`. Supported passes
must preserve the identity of unchanged subtrees so convergence can be detected.

Wrapper extraction can reorder geometric paths. Generate cached inside states
and previous-hit paths for the returned geometry rather than reusing those from
the original representation. Property paths retain their original ordering.
"""
function simplify end

"""
    finalise!(geometry::PhysicalGeometry)

Initialise acceleration structures and derived geometry data after simplification.
"""
function finalise! end

function simplify(geometry::PhysicalGeometry)
    while true
        simplified = simplify(geometry, Val(:extract_common_wrapper))
        if simplified === geometry
            finalise!(geometry)
            return geometry
        end
        geometry = simplified
    end
end

"""
    IntersectionParams{N}(inside, normal, hit_gap)

Store the incident side, surface normal, and gap status of an intersection.
"""
struct IntersectionParams{N}
    inside::Bool
    normal::SVector{N, Float64}
    hit_gap::Bool
end

"""
    size_scale(geometry::PhysicalGeometry)

Return a characteristic geometry length used to choose simulation step sizes.
"""
function size_scale end

"""
    geometry_mesh(geometry; height=nothing, nsamples=100, bounding_box=nothing)

Return renderable surface data, optionally restricting it to `bounding_box`.
`height` and `nsamples` control representations of extended or curved surfaces.
"""
function geometry_mesh end

"""
    distance_to_surface(geometry::PhysicalGeometry, position)

Return the distance from `position` to the geometry's nearest surface.
"""
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

"""
    find_intersection_requires_inside(::Type{<:PhysicalGeometry}) -> Val

Return `Val(true)` when intersection detection needs cached inside-region indices.
"""
find_intersection_requires_inside(::Type{<:PhysicalGeometry}) = Val(false)

"""
    is_liminal(::Type{<:PhysicalGeometry}) -> Bool

Return whether the type represents a liminal geometry with independently placed cells.
"""
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
    contains_mesh(::Type{<:PhysicalGeometry}) -> Bool

Return whether a geometry contains a mesh in its child hierarchy.
"""
function contains_mesh end

"""
    contains_repeat(::Type{<:PhysicalGeometry}) -> Bool

Return whether a geometry's child hierarchy contains a periodic `Repeat`.
"""
function contains_repeat end
contains_repeat(::Type{<:PhysicalGeometry}) = false


"""
    get_intersection_params(geometry::PhysicalGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, indices, isinside=nothing) -> IntersectionParams{N}

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

"""
    get_intersection_params_requires_inside(::Type{<:PhysicalGeometry}) -> Val

Return `Val(true)` when computing intersection parameters needs cached inside indices.
"""
get_intersection_params_requires_inside(::Type{<:PhysicalGeometry}) = Val(false)

"""
    to_child_coordinates(geometry::PhysicalGeometry, object)

Converts the `position`/`bounding_box` from `geometry` coordinates to the coordinates from the child of `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
to_child_coordinates(::PhysicalGeometry{N}, object) where {N} = object

"""
    from_child_coordinates(geometry::PhysicalGeometry, object)

Converts the `position`/`bounding_box` from child of `geometry` coordinates to the coordinates from the `geometry`.

This should be overwritten for transformations. It should not be called for base obstructions.
"""
from_child_coordinates(::PhysicalGeometry{N}, object) where {N} = object

"""
    to_child_coordinates_normal(geometry::PhysicalGeometry, normal)

Convert a surface normal from parent coordinates to child coordinates.
The default leaves the normal unchanged.
"""
to_child_coordinates_normal(::PhysicalGeometry{N}, object) where {N} = object

"""
    from_child_coordinates(geometry::PhysicalGeometry, params::IntersectionParams)

Convert intersection parameters from child to parent coordinates.
"""
from_child_coordinates(::PhysicalGeometry, params::IntersectionParams) = params

"""
    has_inside(::Type{<:PhysicalGeometry}) -> Bool

Return whether a geometry has an inside region.
"""
function has_inside end

"""
    has_single_inside(::Type{<:PhysicalGeometry}) -> Bool

Return whether the geometry represents a single inside region, identified by `()`.
"""
function has_single_inside end

"""
    inside_indices_eltype(::Type{<:PhysicalGeometry})
    inside_indices_eltype(::Type{Union{}})

Return the element type of the array returned by `inside_indices`.
The empty union maps to itself for geometry hierarchies without inside indices.
"""
function inside_indices_eltype end
inside_indices_eltype(::Type{Union{}}) = Union{}

"""
    intersection_type(::Type{<:PhysicalGeometry})
    intersection_type(::Type{Union{}})

Return the surface-index tuple type used by `find_intersection`, excluding its
final incident-side and distance entries.
The empty union maps to itself for geometry hierarchies without surface indices.
"""
function intersection_type end
intersection_type(::Type{Union{}}) = Union{}

"""
    bound_intersection_type(geometry::PhysicalGeometry, density)

Return intersection index types for surfaces with nonzero binding density.
"""
function bound_intersection_type end

"""
    _merge_types(types)

Combine an iterable of types into a union, ignoring empty unions and duplicates.
"""
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

"""
    inside_indices(geometry::PhysicalGeometry, position, intersection=nothing; no_deproject=false)

Return the obstruction indices containing `position`, optionally using a known intersection.
"""
function inside_indices end

"""
    volume_sampling(geometry::PhysicalGeometry, bounding_box, volume_density; no_deproject=false)

Sample positions uniformly throughout `bounding_box`, returning positions and
their cached inside-region indices at the requested number density.
"""
function volume_sampling end

"""
    inside_sampling(geometry::PhysicalGeometry, density; rng=default_rng(), no_deproject=false, bounding_box=nothing)

Sample positions uniformly within a geometry at the requested number density,
returning positions and their cached inside-region indices. When supplied,
`bounding_box` restricts the returned positions.
"""
function inside_sampling end

"""
    inside_sampling_per_obstruction(geometry, density, rng; no_deproject=false, bounding_box=nothing)

Generate inside-volume sampling proposals for one obstruction or a geometry group.
Proposals are not guaranteed to lie within the optional bounding box; callers
must filter them when that restriction is required.
"""
function inside_sampling_per_obstruction end

"""
    projected_surface_area(surface_distribution, direction)

Return the surface area projected perpendicular to `direction`.
"""
function projected_surface_area end

"""
    normal_second_moment(geometry_or_surface_distribution, direction)

Return the direction-dependent second moment of the surface-normal distribution.
"""
function normal_second_moment end

"""
    tortuosity_tensor(liminal_geometry, inside)

Return the diffusion tortuosity tensor for the supplied cached inside state for a liminal geometry.
"""
function tortuosity_tensor end

"""
    inverse_mean_free_path(gliminal_eometry::FixedLiminalGeometry, direction)

Return the inverse mean distance between surface encounters along `direction` for a liminal geometry.
"""
function inverse_mean_free_path end

"""
    _volume_inside_indices(geometry::PhysicalGeometry, position; no_deproject=false)

Determine cached inside-region indices for a sampled position, including the
special cases of geometries with no inside or a single inside region.
"""
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

function inside_sampling_per_obstruction(
    geometry::PhysicalGeometry{N},
    density::Real,
    rng;
    no_deproject=false,
    bounding_box=nothing,
) where {N}
    has_inside(typeof(geometry)) || return SVector{N, Float64}[]
    throw(ArgumentError(
        "inside_sampling_per_obstruction is not implemented for $(typeof(geometry))",
    ))
end

"""
    inside_sampling(geometry::PhysicalGeometry, density; rng=default_rng(), no_deproject=false, bounding_box=nothing)

Sample points uniformly within a geometry, returning positions and inside indices.

When `bounding_box` is provided, retain only samples within that box. The box is
also passed to proposal samplers that need it to bound otherwise infinite
geometries, such as `Repeat`.
"""
function inside_sampling(
    geometry::PhysicalGeometry{N},
    density::Real;
    rng=default_rng(),
    no_deproject=false,
    bounding_box=nothing,
) where {N}
    isfinite(density) && density >= 0 ||
        throw(ArgumentError("inside-sampling density must be finite and nonnegative"))
    if bounding_box isa BoundingBox
        bounding_box = InternalBoundingBox(bounding_box)
    end
    if !isnothing(bounding_box)
        bounding_box isa InternalBoundingBox{N} || throw(ArgumentError(
            "inside-sampling bounding box must be an InternalBoundingBox{$N} or a compatible BoundingBox",
        ))
    end

    proposals = inside_sampling_per_obstruction(
        geometry, density, rng; no_deproject, bounding_box,
    )
    if !isnothing(bounding_box)
        filter!(position -> InternalBoundingBoxes.isinside(bounding_box, position), proposals)
    end
    inside_type = inside_indices_eltype(typeof(geometry))
    positions = SVector{N, Float64}[]
    indices = Vector{inside_type}[]
    for position in proposals
        inside = _volume_inside_indices(geometry, position; no_deproject)
        isempty(inside) && continue
        rand(rng) < inv(length(inside)) || continue
        push!(positions, position)
        push!(indices, inside)
    end
    positions, indices
end

"""
    isinside_single(geometry::PhysicalGeometry, position, previous_intersection=nothing; no_deproject=false)

Return whether `position` lies in the geometry's single inside region.
"""
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

"""
    to_inside_index(geometry::PhysicalGeometry, indices)

Convert a surface intersection path into an inside-region path, or return
`nothing` for geometries without an inside region.
"""
function to_inside_index(geometry::PhysicalGeometry, indices)
    has_inside(typeof(geometry)) || return nothing
    has_single_inside(typeof(geometry)) && return ()
    child, child_indices = get_child(geometry, indices)
    converted = to_inside_index(child, child_indices)
    isnothing(converted) && return nothing
    n_consumed = length(indices) - length(child_indices)
    (indices[1:n_consumed]..., converted...)
end

"""
    random_surface_positions(geometry, density, bounding_box, scale_density; include_gap=false, kwargs...)

Sample surface positions together with initialized collision states, scaling
surface number densities by `scale_density`.
"""
function random_surface_positions end

"""
    InternalBoundingBox(geometry::PhysicalGeometry; kwargs...)

Return the axis-aligned bounding box containing a physical geometry.
"""
function InternalBoundingBox end

"""
    number_of_leaves(geometry::PhysicalGeometry) -> Int

Return the number of leaf surfaces represented by a physical geometry.
"""
function number_of_leaves end

"""
    filter_bounding_box(geometry::PhysicalGeometry, box::InternalBoundingBox, preserve_type::Val=Val(false))

Return the subset of a geometry overlapping `box`, or `nothing`. With `Val(true)`,
preserve child structure during materialization to keep vector child types consistent.
"""
filter_bounding_box(
    geometry::PhysicalGeometry{N}, box::InternalBoundingBox{N}, preserve_type::Val=Val(false),
) where {N} = geometry_from_indices(
    geometry, filter_bounding_box_indices(geometry, box), preserve_type,
)

"""
    filter_bounding_box_indices(geometry, box[, parent_selection])

Select candidates within `box`, optionally refining a parent selection against
the same source geometry. Selection is independent of materialization policy.
Selections preserve original geometry indices and may
contain additional state for inside queries. A selection is always returned;
use `selection_overlaps` to test whether it retains geometry.
"""
function filter_bounding_box_indices end

"""
    geometry_from_indices(geometry, selection, preserve_type)

Materialize a bounding-box selection against its original source geometry.
With `Val(true)`, preserve child structure so vector children have consistent types.
"""
function geometry_from_indices end

"""
    selection_overlaps(selection) -> Bool

Whether a bounding-box selection retains any geometry. Grid construction uses
this predicate to prune empty subtrees without materializing intermediate
geometries. Selection representations with additional state should specialize
this function when they can represent an empty selection without returning
`nothing`.
"""
function selection_overlaps end

selection_overlaps(::Any) = true
selection_overlaps(selection::Bool) = selection

"""
    _bounding_boxes_overlap(first::InternalBoundingBox, second::InternalBoundingBox) -> Bool

Return whether two axis-aligned boxes overlap, including touching boundaries.
"""
function _bounding_boxes_overlap(
    first::InternalBoundingBox{N}, second::InternalBoundingBox{N},
) where {N}
    all(InternalBoundingBoxes.lower(first) .<= InternalBoundingBoxes.upper(second)) &&
        all(InternalBoundingBoxes.lower(second) .<= InternalBoundingBoxes.upper(first))
end

"""
    leaf_bounding_boxes(geometry::PhysicalGeometry; inside=false)

Return conservative bounding boxes for the leaf surfaces of `geometry`, all in
the geometry's coordinate frame. The result mirrors the intersection-index
hierarchy: vector and tuple groups return vectors and tuples of child results,
while transformations and transparent wrappers preserve their child's structure.
Empty groups return empty collections. With `inside=true`, return bounds for
inside regions instead of collision leaves; geometries without an inside return
an empty bounding-box vector while retaining their position in the hierarchy.

A base obstruction returns its `InternalBoundingBox` directly. `FiniteCylinder`
returns a tuple for its side, first cap, and second cap, in surface-index order;
for inside bounds it returns its single full-geometry bounding box. Meshes return
a vector of triangle boxes, including gap-closing triangles, or their full bounding
box for inside queries.

Periodic geometries and dimension-reducing transformations have no finite leaf
bounding boxes in their outer frame. Query their finite child geometry instead.
Likewise, query individual cell geometries for liminal geometry.
"""
function leaf_bounding_boxes end

"""
    _geometry_mesh(geometry; height=nothing, nsamples=100, bounding_box=nothing)

Implement native surface rendering for `geometry` after public bounding boxes
have been converted to the internal representation.
"""
function _geometry_mesh end

"""
    _mesh_result(vertices, triangles)

Package vertices and triangle indices in the named-tuple surface representation.
"""
_mesh_result(vertices, triangles) = (vertices=vertices, triangles=triangles)

"""
    _translate_native(value, displacement)

Translate native surface data, recursively handling collections of surfaces.
"""
_translate_native(value, displacement) = value isa Real ? value + displacement[1] :
    value isa SVector ? value + displacement :
    value isa NamedTuple ? _mesh_result(value.vertices .+ Ref(displacement), value.triangles) :
    [_translate_native(item, displacement) for item in value]

function geometry_mesh(geometry; height=nothing, nsamples=100, bounding_box=nothing)
    bounding_box = bounding_box isa BoundingBox ? InternalBoundingBox(bounding_box) : bounding_box
    _geometry_mesh(geometry; height, nsamples, bounding_box)
end

"""
    show(io::IO, geometry::PhysicalGeometry)

Write a compact representation consisting of the geometry type and an ellipsis.
"""
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
include("liminal_geometries.jl")

import .MonteCarloMeasures: SurfaceEstimate, estimate_volume, estimate_surface
import .Transparents: Transparent, IgnoreOverlapping
import .LiminalGeometries: FixedLiminalGeometry

end
