"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, find_intersection_requires_inside, get_child, get_intersection_params_requires_inside, has_inside, has_single_inside, inside_indices_eltype, intersection_type, bound_intersection_type, _merge_types, isinside_single, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, distance_to_surface, _geometry_mesh
import ..PhysicalGeometries: inside_sampling_per_obstruction
import ..PhysicalGeometries: leaf_bounding_boxes
import ..PhysicalGeometries: filter_bounding_box
import ..PhysicalGeometries: number_of_leaves
import ..PhysicalGeometries: contains_mesh
import ..PhysicalGeometries: contains_repeat
import ..PhysicalGeometries: filter_bounding_box_indices, geometry_from_indices
import ..PhysicalGeometries: selection_overlaps
import ..PhysicalGeometries: has_equivalent_wrapper
import ..PhysicalGeometries: rewrap, rewrap_extracted_wrapper, child_geometry, finalise!
import ..PhysicalGeometries: simplify
import ...Properties: GeometryProperties, GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties
import ...Properties: all_property_values
import ...InsideViews: child_view
import ....BoundingBoxes: BoundingBoxNotSupported

abstract type GroupGeometry{N, P} <: PhysicalGeometry{N} end

contains_mesh(::Type{<:GroupGeometry{N, P}}) where {N, P} = contains_mesh(P)
contains_repeat(::Type{<:GroupGeometry{N, P}}) where {N, P} = contains_repeat(P)

struct GeometryVector{N, P<:PhysicalGeometry{N}, I<:Union{Nothing, Vector{Int}}} <: GroupGeometry{N, P}
    geometries::Vector{P}
    indices::I

    function GeometryVector{N, P, I}(geometries::Vector{P}, indices::I) where {
        N, P<:PhysicalGeometry{N}, I<:Union{Nothing, Vector{Int}},
    }
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        !isnothing(indices) && length(indices) != length(geometries) && throw(ArgumentError(
            "indices and geometries must have the same length",
        ))
        !isnothing(indices) && !allunique(indices) && throw(ArgumentError(
            "group indices must be unique",
        ))
        new{N, P, I}(geometries, indices)
    end
end

number_of_leaves(group::GeometryVector) =
    sum(number_of_leaves, group_geometries(group); init=0)

function indices end

child_geometry(group::GroupGeometry) = throw(ArgumentError(
    "geometry group $(typeof(group)) has multiple children; use group_geometries instead",
))

function inside_sampling_per_obstruction(
    group::GroupGeometry{N}, density::Real, rng; no_deproject=false, bounding_box=nothing,
) where {N}
    proposals = SVector{N, Float64}[]
    for child in group_geometries(group; include_gap=false)
        append!(proposals, inside_sampling_per_obstruction(
            child, density, rng; no_deproject, bounding_box,
        ))
    end
    proposals
end

get_intersection_params_requires_inside(::Type{<:GroupGeometry{N, P}}) where {N, P} =
    get_intersection_params_requires_inside(P)

_trait_as_val(::Val{value}) where {value} = Val(value)
_trait_as_val(value::Bool) = Val(value)

_tuple_trait_any(::F, ::Type{Tuple{}}) where {F} = Val(false)

_tuple_trait_any(trait, ::Type{P}) where {P<:Tuple} =
    _tuple_trait_any(
        trait,
        _trait_as_val(trait(Base.tuple_type_head(P))),
        Base.tuple_type_tail(P),
    )

_tuple_trait_any(::F, ::Val{true}, ::Type) where {F} = Val(true)
_tuple_trait_any(::F, ::Val{false}, ::Type{Tuple{}}) where {F} = Val(false)

_tuple_trait_any(trait, ::Val{false}, ::Type{P}) where {P<:Tuple} =
    _tuple_trait_any(
        trait,
        _trait_as_val(trait(Base.tuple_type_head(P))),
        Base.tuple_type_tail(P),
    )

function Base.show(io::IO, ::Type{T}) where {T <: GeometryVector}
    print(io, "GeometryVector")
    T isa UnionAll && return
    length(T.parameters) < 2 && return
    print(io, "{")
    show(io, T.parameters[2])
    print(io, "}")
end

group_geometries(group::GeometryVector; include_gap=true) = group.geometries
indices(group::GeometryVector) = isnothing(group.indices) ? eachindex(group.geometries) : group.indices

function _indexed_children(group::GroupGeometry; include_gap=true)
    zip(indices(group), group_geometries(group; include_gap))
end

struct GroupSelection{G, S, B}
    original::Ref{G}
    child_selections::Vector{S}
    box::B
end

selection_overlaps(selection::GroupSelection) = any(selection_overlaps, selection.child_selections)
selection_overlaps(selection::GroupSelection{G, Bool}) where {G} = iszero(sum(selection.child_selections))

function filter_bounding_box_indices(group::GroupGeometry, box)
    children = [filter_bounding_box_indices(child, box) for child in group_geometries(group)]
    GroupSelection(Ref(group), children, box)
end

function filter_bounding_box_indices(group::GroupGeometry, box, selection::GroupSelection{G, S}) where {G, S}
    children = Vector{S}(undef, length(selection.child_selections))
    group_members = group_geometries(group)
    for index in eachindex(children)
        @inbounds children[index] = filter_bounding_box_indices(group_members[index], box, selection.child_selections[index])
    end
    GroupSelection(selection.original, children, box)
end

function filter_bounding_box_indices(group::GroupGeometry, box, selection::GroupSelection{G, Bool}) where {G}
    if iszero(sum(selection.child_selections))
        return selection
    end
    children = Vector{Bool}(undef, length(selection.child_selections))
    group_members = group_geometries(group)
    for index in eachindex(group_members)
        @inbounds children[index] = selection.child_selections[index] && filter_bounding_box_indices(group_members[index], box, selection.child_selections[index])
    end
    GroupSelection(selection.original, children, box)
end


function geometry_from_indices(
    group::GeometryVector, selection::GroupSelection, ::Val{preserve_type},
) where {preserve_type}
    retained = findall(selection_overlaps, selection.child_selections)
    original = selection.original[]
    selected_indices = [indices(original)[position] for position in retained]
    isempty(selected_indices) && !preserve_type && return nothing
    if isempty(selected_indices)
        source_children = group_geometries(original)
        isempty(source_children) && return rewrap(group, eltype(source_children)[]; indices=Int[])
        first_child = first(source_children)
        child_selection = filter_bounding_box_indices(first_child, selection.box)
        child_type = typeof(geometry_from_indices(first_child, child_selection, Val(true)))
        return rewrap(group, child_type[]; indices=Int[])
    end
    selected_selections = selection.child_selections[retained]
    children = map(retained, selected_selections) do position, child_selection
        geometry_from_indices(group_geometries(original)[position], child_selection, Val(preserve_type))
    end
    Child = typeof(first(children))
    if all(child -> typeof(child) === Child, children)
        typed_children = children isa Vector{Child} ? children : Vector{Child}(children)
        return rewrap(group, typed_children; indices=selected_indices)
    end
    preserve_type && throw(ArgumentError(
        "type-preserving filtering produced an inconsistent GeometryVector child type",
    ))
    geometry_from_indices(group, selection, Val(true))
end

function distance_to_surface(group::GroupGeometry{N}, position::SVector{N, Float64}) where {N}
    minimum(
        (distance_to_surface(child, position)
         for child in group_geometries(group; include_gap=false));
        init=Inf,
    )
end

child_type(::Type{<:GroupGeometry{N, P}}) where {N, P} = P

function _prepend_type(::Type{Prefix}, ::Type{T}) where {Prefix, T}
    T === Union{} && return Union{}
    T isa Union && return Union{
        (_prepend_type(Prefix, element) for element in Base.uniontypes(T))...,
    }
    T <: Tuple && return Tuple{Prefix, T.parameters...}
    throw(MethodError(_prepend_type, (Type{Prefix}, Type{T})))
end

function _append_type(::Type{T}, ::Type{Suffix}) where {T, Suffix}
    T === Union{} && return Union{}
    T isa Union && return Union{
        (_append_type(element, Suffix) for element in Base.uniontypes(T))...,
    }
    T <: Tuple || throw(MethodError(_append_type, (Type{T}, Type{Suffix})))
    Tuple{T.parameters..., Suffix}
end

_union_val(::Val{A}, ::Val{B}) where {A, B} = Val(Union{A, B})
_type_from_val(::Val{T}) where {T} = T

_inside_indices_union(::Type{Union{A, B}}) where {A, B} = _union_val(
    Val(inside_indices_eltype(A)),
    Val(inside_indices_eltype(B)),
)

_intersection_type_union(::Type{Union{A, B}}) where {A, B} = _union_val(
    Val(intersection_type(A)),
    Val(intersection_type(B)),
)

_inside_indices_tuple(::Type{Tuple{}}) = Val(Union{})
_inside_indices_tuple(::Type{P}) where {P<:Tuple} = _union_val(
    Val(_prepend_type(Int, inside_indices_eltype(Base.tuple_type_head(P)))),
    _inside_indices_tuple(Base.tuple_type_tail(P)),
)

_intersection_type_tuple(::Type{Tuple{}}) = Val(Union{})
_intersection_type_tuple(::Type{P}) where {P<:Tuple} = _union_val(
    Val(_prepend_type(Int, intersection_type(Base.tuple_type_head(P)))),
    _intersection_type_tuple(Base.tuple_type_tail(P)),
)

inside_indices_eltype(::Type{Union{A, B}}) where {A, B} =
    _type_from_val(_inside_indices_union(Union{A, B}))

intersection_type(::Type{Union{A, B}}) where {A, B} =
    _type_from_val(_intersection_type_union(Union{A, B}))

inside_indices_eltype(::Type{<:GroupGeometry{N, P}}) where {N, P} =
    _prepend_type(Int, inside_indices_eltype(P))

intersection_type(::Type{<:GroupGeometry{N, P}}) where {N, P} =
    _prepend_type(Int, intersection_type(P))

function bound_intersection_type(geometry::GroupGeometry, density)
    _merge_types(
        (
            _prepend_type(
                Int,
                bound_intersection_type(child, child_density),
            )
            for (index, child) in _indexed_children(geometry)
            for child_density in (_density_child(density, index),)
            if any(!iszero, all_property_values(child_density))
        )
    )
end


find_intersection_requires_inside(::Type{<:GroupGeometry{N, P}}) where {N, P} =
    find_intersection_requires_inside(P)

function find_intersection(group::GroupGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, previous_hit=nothing, inside=nothing, candidate::Union{Nothing, Tuple}=nothing) where {N}
    if !isnothing(candidate)
        isempty(candidate) && throw(ArgumentError("a group candidate must identify a child"))
        index = first(candidate)
        child = group_geometries(group)[_stored_index(group, index)]
        child_previous = isnothing(previous_hit) || previous_hit[1] != index ? nothing : Base.tail(previous_hit)
        child_inside = child_view(find_intersection_requires_inside(typeof(child)), inside, index)
        hit = find_intersection(child, start, dest, child_previous, child_inside, Base.tail(candidate))
        return isnothing(hit) ? nothing : (index, hit...)
    end
    current = nothing
    for (index, candidate, dist_all_checked) in intersection_candidates(group, start, dest)
        if !isnothing(current) && current[end] < dist_all_checked
            return current
        end
        candidate_previous_hit = isnothing(previous_hit) || previous_hit[1] != index ? nothing : previous_hit[2:end]
        candidate_inside = child_view(
            find_intersection_requires_inside(typeof(candidate)), inside, index
        )
        intersect = find_intersection(candidate, start, dest, candidate_previous_hit, candidate_inside)
        if isnothing(intersect)
            continue
        end
        if isnothing(current) || intersect[end] < current[end]
            current = (index, intersect...)
        end
    end
    return current
end

_is_fully_concrete(::Type{T}) where {T} =
    isconcretetype(T) && all(_child_geometry_types_are_concrete, T.parameters)

function _child_geometry_types_are_concrete(::Type{T}) where {T}
    T isa Union && return all(_child_geometry_types_are_concrete, Base.uniontypes(T))
    T <: PhysicalGeometry && return _is_fully_concrete(T)
    T <: Tuple && return all(_child_geometry_types_are_concrete, T.parameters)
    true
end
_child_geometry_types_are_concrete(::Any) = true

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}, I<:Union{Nothing, Vector{Int}}} <: GroupGeometry{N, P}
    geometries::P
    indices::I

    function GeometryTuple{N, P, I}(geometries::P, indices::I) where {
        N, P<:Tuple{Vararg{PhysicalGeometry{N}}}, I<:Union{Nothing, Vector{Int}},
    }
        !isnothing(indices) && length(indices) != length(geometries) && throw(ArgumentError(
            "indices and geometries must have the same length",
        ))
        !isnothing(indices) && !allunique(indices) && throw(ArgumentError(
            "group indices must be unique",
        ))
        new{N, P, I}(geometries, indices)
    end
end

contains_mesh(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _tuple_trait_any(contains_mesh, P) === Val(true)
contains_repeat(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    any(contains_repeat, P.parameters)

function geometry_from_indices(
    group::GeometryTuple, selection::GroupSelection, ::Val{preserve_type},
) where {preserve_type}
    preserve_type || begin
        retained = findall(selection_overlaps, selection.child_selections)
        original = selection.original[]
        children = map(retained) do position
            geometry_from_indices(
                group_geometries(original)[position],
                selection.child_selections[position],
                Val(false),
            )
        end
        isempty(children) && return nothing
        selected_indices = [indices(original)[position] for position in retained]
        return rewrap(group, Tuple(children); indices=selected_indices)
    end
    original = selection.original[]
    children = map(enumerate(original.geometries)) do (position, child)
        child_selection = selection.child_selections[position]
        selection_overlaps(child_selection) ?
            geometry_from_indices(child, child_selection, Val(true)) :
            geometry_from_indices(child, filter_bounding_box_indices(child, selection.box), Val(true))
    end
    child_tuple = Tuple(children)
    all(child_tuple[i] === original.geometries[i] for i in eachindex(child_tuple)) && return original
    rewrap(original, child_tuple; indices=original.indices)
end

number_of_leaves(group::GeometryTuple) =
    sum(number_of_leaves, group.geometries; init=0)

function find_intersection_requires_inside(::Type{<:GeometryTuple{N, P}}) where {N, P}
    _tuple_trait_any(find_intersection_requires_inside, P)
end

get_intersection_params_requires_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _tuple_trait_any(get_intersection_params_requires_inside, P)

group_geometries(group::GeometryTuple; include_gap=true) = group.geometries
indices(group::GeometryTuple) = isnothing(group.indices) ? Base.OneTo(length(group.geometries)) : group.indices

function _stored_index(group::GroupGeometry, index::Int)
    group_indices = indices(group)
    if group_indices isa Base.OneTo
        checkbounds(Bool, group_indices, index) || throw(BoundsError(group, index))
        return index
    end
    position = findfirst(==(index), group_indices)
    isnothing(position) && throw(BoundsError(group, index))
    position
end

_indexed_children(group::GeometryTuple; include_gap=true) = zip(indices(group), group.geometries)

function has_equivalent_wrapper(group::GroupGeometry)
    children = group_geometries(group)
    isempty(children) && return false
    reference = first(children)
    all(child -> has_equivalent_wrapper(reference, child), children)
end

rewrap(group::GeometryVector, children::AbstractVector{<:PhysicalGeometry}; indices=group.indices) =
    GeometryVector(children; indices)
rewrap(group::GeometryTuple, children::Tuple; indices=group.indices) =
    GeometryTuple(children; indices)

function finalise!(group::Union{GeometryVector, GeometryTuple})
    foreach(finalise!, group_geometries(group))
    nothing
end

# Return the rebuilt node, or the original node when child changes are rejected.
function _simplify_children(group::GeometryVector, pass::Val)
    isempty(group) && return group
    children = map(child -> simplify(child, pass), group.geometries)
    child_type = typeof(first(children))
    _is_fully_concrete(child_type) && all(child -> typeof(child) === child_type, children) || return group
    unchanged = all(index -> children[index] === group.geometries[index], eachindex(children))
    unchanged ? group : rewrap(group, child_type[children...])
end

function _simplify_children(group::GeometryTuple, pass::Val)
    children = map(child -> simplify(child, pass), group.geometries)
    unchanged = all(index -> children[index] === group.geometries[index], eachindex(children))
    unchanged ? group : rewrap(group, children)
end

simplify(group::Union{GeometryVector, GeometryTuple}, pass::Val) = _simplify_children(group, pass)

function simplify(group::Union{GeometryVector, GeometryTuple}, pass::Val{:extract_common_wrapper})
    rebuilt = _simplify_children(group, pass)
    has_equivalent_wrapper(rebuilt) || return rebuilt
    children = map(child_geometry, rebuilt.geometries)
    rewrap_extracted_wrapper(first(rebuilt.geometries), rebuilt.geometries, rewrap(rebuilt, children))
end

leaf_bounding_boxes(group::GroupGeometry; inside=false) = map(
    child -> leaf_bounding_boxes(child; inside), group_geometries(group),
)

function get_child(group::GroupGeometry, child_indices::Tuple)
    child_index = child_indices[1]
    group_geometries(group)[_stored_index(group, child_index)], Base.tail(child_indices)
end

child_type(::Type{<:GeometryTuple{N, P}}) where {N, P} = Union{P.parameters...}

inside_indices_eltype(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _type_from_val(_inside_indices_tuple(P))

intersection_type(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _type_from_val(_intersection_type_tuple(P))

function bound_intersection_type(geometry::GeometryTuple, density)
    _merge_types(
        (
            _prepend_type(
                Int,
                bound_intersection_type(child, _density_child(density, index)),
            )
            for (index, child) in _indexed_children(geometry)
            if any(!iszero, all_property_values(_density_child(density, index)))
        )
    )
end

function inside_candidates end

inside_candidates(group::GeometryVector, position) = _indexed_children(group)

inside_candidates(group::GeometryTuple, position) = _indexed_children(group)

function isinside_single(
    group::GroupGeometry{N},
    position::SVector{N, Float64},
    previous_intersection=nothing;
    candidate=nothing,
    no_deproject=false,
    kwargs...,
) where {N}
    isnothing(candidate) && throw(ArgumentError(
        "isinside_single for a geometry group requires a candidate index",
    ))
    candidate isa Tuple && !isempty(candidate) || throw(ArgumentError(
        "a geometry-group candidate must be a non-empty index tuple",
    ))

    child_index = candidate[1]
    child, child_candidate = get_child(group, candidate)
    has_inside(typeof(child)) || throw(ArgumentError(
        "candidate index selects a geometry without an inside region",
    ))

    child_intersection = !isnothing(previous_intersection) &&
        previous_intersection[1] == child_index ? previous_intersection[2:end] : nothing
    if has_single_inside(typeof(child))
        isempty(child_candidate) || throw(ArgumentError(
            "candidate index contains extra indices for a single-inside geometry",
        ))
        return isinside_single(child, position, child_intersection; no_deproject, kwargs...)
    end
    isinside_single(
        child, position, child_intersection;
        candidate=child_candidate, no_deproject, kwargs...,
    )
end

intersection_candidates(group::GroupGeometry, start, destination) =
    ((index, child, 0.0) for (index, child) in _indexed_children(group))

has_inside(::Type{<:GroupGeometry{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _tuple_trait_any(has_inside, P) === Val(true)
has_single_inside(::Type{<:GroupGeometry}) = false

function InternalBoundingBox(geometry::GroupGeometry{N, P}; kwargs...) where {N, P}
    isempty(geometry) && throw(BoundingBoxNotSupported(
        "cannot construct a bounding box for an empty geometry group",
    ))
    boxes = [InternalBoundingBoxes.InternalBoundingBox(child; kwargs...) for child in geometry]
    lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in boxes))
    upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in boxes))
    InternalBoundingBoxes.InternalBoundingBox(
        lower_bound,
        upper_bound,
    )
end

function inside_indices_for_any_type(geometry, position, intersection; no_deproject=false)
    geometry_type = typeof(geometry)
    has_inside(geometry_type) || return Tuple{}[]
    if has_single_inside(geometry_type)
        return isinside_single(geometry, position, intersection; no_deproject) ?
            [()] : Tuple{}[]
    end
    return inside_indices(geometry, position, intersection; no_deproject)
end

function inside_indices(
    geometry::GroupGeometry{N, P},
    position::SVector{N, Float64},
    intersection=nothing,
    ; no_deproject=false,
) where {N, P}
    indices = Vector{inside_indices_eltype(typeof(geometry))}()
    if !has_inside(P)
        return indices
    end
    single_inside = has_single_inside(P)
    for (child_index, child) in inside_candidates(geometry, position)
        child_intersection = !isnothing(intersection) && intersection[1] == child_index ? intersection[2:end] : nothing
        if single_inside
            if isinside_single(child, position, child_intersection; no_deproject)
                push!(indices, (child_index, ))
            end
        else
            append!(indices, [(child_index, new_indices...) for new_indices in inside_indices(
                child, position, child_intersection; no_deproject,
            )])
        end
    end
    return indices
end

function inside_indices(
    geometry::GeometryTuple{N},
    position::SVector{N, Float64},
    intersection=nothing;
    no_deproject=false,
) where {N}
    indices = Vector{inside_indices_eltype(typeof(geometry))}()
    for (child_index, child) in _indexed_children(geometry)
        child_intersection = !isnothing(intersection) && intersection[1] == child_index ? intersection[2:end] : nothing
        append!(indices, [(child_index, new_indices...) for new_indices in inside_indices_for_any_type(
            child, position, child_intersection; no_deproject,
        )])
    end
    indices
end

GeometryVector{N, P}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P, Nothing}(geometries, nothing)

GeometryVector{N, P}(geometries::Vector{P}, indices::Vector{Int}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P, Vector{Int}}(geometries, indices)

GeometryVector{N}(geometries::Vector{P}; indices=nothing) where {N, P<:PhysicalGeometry{N}} =
    isnothing(indices) ? GeometryVector{N, P}(geometries) : GeometryVector{N, P}(geometries, indices)

GeometryVector(geometries::Vector{P}; indices=nothing) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N}(geometries; indices)

GeometryTuple{N, P}(geometries::P) where {N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} =
    GeometryTuple{N, P, Nothing}(geometries, nothing)

GeometryTuple{N, P}(geometries::P, indices::Vector{Int}) where {N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} =
    GeometryTuple{N, P, Vector{Int}}(geometries, indices)

GeometryTuple{N}(geometries::Tuple{Vararg{PhysicalGeometry{N}}}; indices=nothing) where {N} =
    isnothing(indices) ? GeometryTuple{N, typeof(geometries)}(geometries) :
        GeometryTuple{N, typeof(geometries)}(geometries, Vector{Int}(indices))

GeometryTuple(geometries::Tuple{Vararg{PhysicalGeometry{N}}}; indices=nothing) where {N} =
    GeometryTuple{N}(geometries; indices)

Base.size(geometry::GroupGeometry) = size(group_geometries(geometry))
Base.axes(geometry::GroupGeometry) = axes(group_geometries(geometry))
Base.length(geometry::GroupGeometry) = length(group_geometries(geometry))
Base.firstindex(geometry::GroupGeometry) = firstindex(group_geometries(geometry))
Base.lastindex(geometry::GroupGeometry) = lastindex(group_geometries(geometry))
Base.getindex(geometry::GroupGeometry, index...) = getindex(group_geometries(geometry), index...)
Base.iterate(geometry::GroupGeometry, state...) = iterate(group_geometries(geometry), state...)
Base.eltype(::Type{<:GroupGeometry{N, P}}) where {N, P} = P
Base.setindex!(geometry::GeometryVector, value, index...) = setindex!(group_geometries(geometry), value, index...)
Base.eltype(::Type{GeometryTuple{N, P}}) where {N, P} = eltype(P)
Base.Tuple(geometry::GeometryTuple) = group_geometries(geometry)

function Base.show(io::IO, ::Type{T}) where {T <: GeometryTuple}
    print(io, "GeometryTuple")
    T isa UnionAll && return
    print(io, "{")
    for (index, child) in enumerate(T.parameters[2].parameters)
        index > 1 && print(io, ", ")
        show(io, child)
    end
    print(io, "}")
end

_density_child(density::GeometryLeafProperties, ::Int) = density
_density_child(density::GeometryVectorProperties, index::Int) = density.properties[index]
_density_child(density::GeometryTupleProperties, index::Int) = density.properties[index]
_density_child(density::GeometryProperties, ::SVector) = density

function _combine(draws, ::Val{N}, ::Type{IndexType}) where {N, IndexType}
    draws = collect(draws)
    npositions = sum(length(draw[1]) for draw in draws; init=0)
    nindices = sum(length(draw[2]) for draw in draws; init=0)
    positions = Vector{SVector{N, Float64}}(undef, npositions)
    indices = Vector{IndexType}(undef, nindices)

    position_start = 1
    index_start = 1
    for (child_positions, child_indices) in draws
        position_end = position_start + length(child_positions) - 1
        index_end = index_start + length(child_indices) - 1
        copyto!(positions, position_start, child_positions, firstindex(child_positions), length(child_positions))
        copyto!(indices, index_start, child_indices, firstindex(child_indices), length(child_indices))
        position_start = position_end + 1
        index_start = index_end + 1
    end
    positions, indices
end

function random_surface_positions(geometry::GroupGeometry{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density;
    include_gap=false, kwargs...
) where {N}
    IndexType = _append_type(intersection_type(typeof(geometry)), Bool)
    (iszero(scale_density) ||
     (density isa GeometryLeafProperties && iszero(density.value))) &&
        return SVector{N, Float64}[], IndexType[]

    _combine((let
        values = random_surface_positions(
            child,
            _density_child(density, index),
            bounding_box,
            scale_density;
            include_gap, kwargs...,
        )
        (values[1], [(index, child_index...) for child_index in values[2]])
     end for (index, child) in _indexed_children(geometry; include_gap)),
        Val(N),
        IndexType,
    )
end

size_scale(geometry::GroupGeometry) = isempty(geometry) ? Inf : minimum(size_scale, geometry)

_geometry_mesh(geometry::GroupGeometry; kwargs...) =
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])

end
