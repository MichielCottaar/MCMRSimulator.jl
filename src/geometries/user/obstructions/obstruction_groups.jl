"""
Defines `ObstructionGroupData`, `ObstructionGroup`, and `IndexedObstruction`.
"""
module ObstructionGroups

import ..ObstructionTypes: ObstructionType, key_value_pairs
import ..Fields: Field, FieldValue, isglobal, value_as_vector
import ....Internal: all_geometry_ids

struct ObstructionGroupData{N}
    type::ObstructionType{N}
    field_values::Dict
    unique_keys::Vector{Symbol}
    n_obstructions::Int
end

"""A group or indexed view over shared obstruction data."""
struct ObstructionGroup{N, I <: AbstractVector{Int}}
    data::ObstructionGroupData{N}
    indices::I
end

struct IndexedObstruction{N}
    data::ObstructionGroupData{N}
    index::Int
end

struct FieldValueView{F <: FieldValue, I <: AbstractVector{Int}}
    field_value::F
    indices::I
    n_data::Int
end

_full_group(group::ObstructionGroup) = group.indices == Base.OneTo(group.data.n_obstructions) ||
    group.indices == (1:group.data.n_obstructions)

geometry_id(group::ObstructionGroup) = objectid(group.data.field_values)
geometry_id(obstruction::IndexedObstruction) = (objectid(obstruction.data.field_values), obstruction.index)

const ObstructionGeometryID = Tuple{UInt, Int}

all_geometry_ids(obstruction::IndexedObstruction) = Set{ObstructionGeometryID}((geometry_id(obstruction),))
all_geometry_ids(group::ObstructionGroup) = Set{ObstructionGeometryID}(
    geometry_id(group[index]) for index in 1:length(group)
)

function _union_geometry_ids(geometries)
    ids = Set{ObstructionGeometryID}()
    for geometry in geometries
        union!(ids, all_geometry_ids(geometry))
    end
    ids
end

all_geometry_ids(geometries::AbstractArray) = _union_geometry_ids(geometries)
all_geometry_ids(geometries::Tuple) = _union_geometry_ids(geometries)

Base.length(group::ObstructionGroup) = length(group.indices)
Base.getindex(group::ObstructionGroup, ::Colon) = group
Base.getindex(group::ObstructionGroup, index::Int) = IndexedObstruction(group.data, group.indices[index])
IndexedObstruction(group::ObstructionGroup, index::Int) = group[index]
Base.getindex(group::ObstructionGroup, indices::AbstractVector{<:Integer}) =
    ObstructionGroup(group.data, Int[group.indices[index] for index in indices])
Base.getindex(group::ObstructionGroup, mask::AbstractVector{Bool}) = group[findall(mask)]
Base.getindex(group::ObstructionGroup, indices::AbstractRange{<:Integer}) =
    ObstructionGroup(group.data, Int[group.indices[index] for index in indices])

Base.iterate(group::ObstructionGroup) = iterate(group, 1)
function Base.iterate(group::ObstructionGroup, state::Int)
    state > length(group) && return nothing
    group[state], state + 1
end

Base.propertynames(group::ObstructionGroup) = group.data.unique_keys
Base.propertynames(obstruction::IndexedObstruction) = obstruction.data.unique_keys

function Base.getproperty(group::ObstructionGroup, name::Symbol)
    name === :data && return getfield(group, :data)
    name === :indices && return getfield(group, :indices)
    data = getfield(group, :data)
    name === :n_obstructions && return length(group)
    name in (:type, :field_values, :unique_keys) && return getproperty(data, name)
    if haskey(data.field_values, name)
        field_value = data.field_values[name]
        return _full_group(group) ? field_value : FieldValueView(field_value, group.indices, data.n_obstructions)
    end
    getproperty(data, name)
end

function Base.setproperty!(group::ObstructionGroup, name::Symbol, value)
    data = getfield(group, :data)
    haskey(data.field_values, name) || return setproperty!(data, name, value)
    field_value = data.field_values[name]
    if _full_group(group) || field_value.field.only_group
        field_value.value = value
        return value
    end
    values = value isa AbstractVector ? value : fill(value, length(group))
    length(values) == length(group) || throw(ArgumentError(
        "Expected one value for each of the $(length(group)) obstructions in this view, not $(length(values)).",
    ))
    for (index, selected_index) in enumerate(group.indices)
        setindex!(field_value, values[index], selected_index)
    end
    value
end

function Base.getproperty(obstruction::IndexedObstruction, name::Symbol)
    name in (:data, :index) && return getfield(obstruction, name)
    data = getfield(obstruction, :data)
    if haskey(data.field_values, name)
        return data.field_values[name][getfield(obstruction, :index)]
    end
    getproperty(data, name)
end

function Base.setproperty!(obstruction::IndexedObstruction, name::Symbol, value)
    data = getfield(obstruction, :data)
    haskey(data.field_values, name) || return setproperty!(data, name, value)
    setindex!(data.field_values[name], value, getfield(obstruction, :index))
end

function Base.getproperty(view::FieldValueView, name::Symbol)
    field_value = getfield(view, :field_value)
    name === :field_value && return field_value
    name === :indices && return getfield(view, :indices)
    name === :n_data && return getfield(view, :n_data)
    name === :field && return field_value.field
    name === :category && return field_value.category
    name === :n_obstructions && return length(view)
    if name === :value
        if isglobal(field_value) || length(field_value) != view.n_data
            return field_value.value
        end
        return Base.view(field_value.value, view.indices)
    end
    getproperty(field_value, name)
end

function Base.setproperty!(view::FieldValueView, name::Symbol, value)
    name === :value || throw(ArgumentError("Cannot set FieldValueView.$name"))
    field_value = getfield(view, :field_value)
    if field_value.field.only_group
        field_value.value = value
        return value
    end
    values = value isa AbstractVector ? value : fill(value, length(view))
    length(values) == length(view) || throw(ArgumentError(
        "Expected one value for each of the $(length(view)) obstructions in this view, not $(length(values)).",
    ))
    for (index, selected_index) in enumerate(view.indices)
        setindex!(field_value, values[index], selected_index)
    end
    value
end

isglobal(view::FieldValueView) = isglobal(view.field_value)
Base.length(view::FieldValueView) = length(view.indices)
function Base.getindex(view::FieldValueView, index::Int)
    view.field_value[view.indices[index]]
end
Base.getindex(view::FieldValueView, indices::AbstractVector{<:Integer}) =
    FieldValueView(view.field_value, Int[view.indices[index] for index in indices], view.n_data)
Base.getindex(view::FieldValueView, mask::AbstractVector{Bool}) = view[findall(mask)]
Base.getindex(view::FieldValueView, indices::AbstractRange{<:Integer}) =
    FieldValueView(view.field_value, Int[view.indices[index] for index in indices], view.n_data)
Base.getindex(view::FieldValueView, ::Colon) = view
function Base.setindex!(view::FieldValueView, value, index::Int)
    setindex!(view.field_value, value, view.indices[index])
end
function Base.setindex!(view::FieldValueView, values::AbstractVector, indices::Union{AbstractVector{<:Integer}, AbstractRange{<:Integer}})
    length(values) == length(indices) || throw(DimensionMismatch("number of values must match the selected field entries"))
    for (value, index) in zip(values, indices)
        setindex!(view, value, index)
    end
    values
end
Base.setindex!(view::FieldValueView, values::AbstractVector, mask::AbstractVector{Bool}) = setindex!(view, values, findall(mask))
value_as_vector(view::FieldValueView) = isglobal(view) ? fill(view.value, length(view)) : collect(view.value)

"""Materialize a group view while preserving global fields and selected source values."""
function materialize(group::ObstructionGroup)
    _full_group(group) && return group
    data = group.data
    n_selected = length(group)
    field_values = Dict{Symbol, FieldValue}()
    copied_values = IdDict{FieldValue, FieldValue}()
    for (name, field_value) in data.field_values
        if haskey(copied_values, field_value)
            field_values[name] = copied_values[field_value]
            continue
        end
        value = field_value.value
        n_values = field_value.n_obstructions
        if !isglobal(field_value) && n_values == data.n_obstructions
            value = value[group.indices]
            n_values = n_selected
        end
        copied = isnothing(field_value.category) ?
            FieldValue(field_value.field, n_values, value) :
            FieldValue(field_value.field, n_values, field_value.category, value)
        field_values[name] = copied
        copied_values[field_value] = copied
    end
    ObstructionGroup(ObstructionGroupData(data.type, field_values, copy(data.unique_keys), n_selected), collect(1:n_selected))
end

"""A single obstruction from an [`ObstructionGroup`](@ref)."""
function Base.show(io::IO, obstruction::IndexedObstruction)
    data = obstruction.data
    print(io, "$(data.type.singular) at index $(obstruction.index) (")
    already_done = Set()
    for (key, field_value) in data.field_values
        field_value in already_done && continue
        push!(already_done, field_value)
        value = field_value[obstruction.index]
        isnothing(value) || print(io, "$key = $value, ")
    end
    print(io, ")")
end

function Base.show(io::IO, group::ObstructionGroup)
    data = group.data
    print(io, "$(length(group)) $(lowercase(String(data.type.plural)))(")
    not_set = Symbol[]
    required_missing = Symbol[]
    already_done = IdDict{Any, Nothing}()
    for key in data.unique_keys
        field_value = getproperty(group, key)
        source = field_value isa FieldValueView ? field_value.field_value : field_value
        haskey(already_done, source) && continue
        already_done[source] = nothing
        if isnothing(field_value.value)
            if field_value.field.required
                push!(required_missing, key)
            else
                push!(not_set, key)
            end
        elseif (isglobal(field_value) && key != :vertices) || length(group) < 5
            print(io, "$key = $(field_value.value), ")
        else
            print(io, "$key = <unique values>, ")
        end
    end
    print(io, ")")
    if !isempty(required_missing)
        print(io, "\nMissing required fields: $(join(required_missing, ", ", " and "))")
    end
    if !isempty(not_set)
        print(io, "\nNot set: $(join(not_set, ", ", " and "))")
    end
end

"""Return the number of enclosed volumes represented by a group view."""
function nvolumes(group::ObstructionGroup)
    data = group.data
    isempty(data.type.volumes) && error("Getting the number of volumes is not defined for $(data.type)")
    data.type.group_volumes ? (isempty(group) ? 0 : 1) : length(group)
end

function ObstructionGroup(type::ObstructionType; number=nothing, kwargs...)
    (unique_keys, field_values) = key_value_pairs(type, isnothing(number) ? 0 : number)
    for (key, value) in kwargs
        field_values[key].value = value
    end

    if isnothing(number)
        for field_value in values(field_values)
            if field_value.field.per_volume && type.group_volumes
                continue
            end
            if !isglobal(field_value)
                number = length(field_value.value)
                break
            end
        end
        isnothing(number) && (number = 1)
    end

    data = ObstructionGroupData(type, field_values, unique_keys, number)
    group = ObstructionGroup(data, collect(1:number))
    for field_value in values(field_values)
        if field_value.field.per_volume && type.group_volumes
            field_value.n_obstructions = nvolumes(group)
        else
            field_value.n_obstructions = number
        end
        if !isglobal(field_value) && length(field_value.value) != field_value.n_obstructions
            error("Incorrect number of values for $field_value. Expected $(field_value.n_obstructions); got $(length(field_value.value)).")
        end
    end
    group
end

end
