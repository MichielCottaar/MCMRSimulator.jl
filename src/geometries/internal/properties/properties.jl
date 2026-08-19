"""MRI-property storage and lookup for physical geometries."""
module Properties

import StaticArrays: SVector

"""
MRI properties assigned to obstructions or groups of obstructions of type `S`.
"""
abstract type GeometryProperties{S} end

Base.eltype(::GeometryProperties{S}) where {S} = S

struct GeometryLeafProperties{S} <: GeometryProperties{S}
    value :: S
end

struct GeometryVectorProperties{S, G<:GeometryProperties{S}} <: GeometryProperties{S}
    properties::Vector{G}
end

struct GeometryTupleProperties{S, P<:Tuple} <: GeometryProperties{S}
    properties::P
end

function _as_properties(value)
    value isa GeometryProperties && return value
    GeometryLeafProperties(value)
end

function GeometryVectorProperties(properties::AbstractVector)
    children = [_as_properties(property) for property in properties]
    isempty(children) && throw(ArgumentError("cannot infer a property type from an empty vector"))
    G = typeof(children[1])
    S = eltype(children[1])
    isconcretetype(S) || throw(ArgumentError("property type must be concrete"))
    all(typeof(property) === G for property in children) ||
        throw(ArgumentError("all vector children must have the same concrete type"))
    all(eltype(property) === S for property in children) ||
        throw(ArgumentError("all property children must have the same concrete type"))
    GeometryVectorProperties{S, G}(G[children...])
end

function GeometryTupleProperties(properties::Tuple)
    children = map(_as_properties, properties)
    isempty(children) && throw(ArgumentError("cannot infer a property type from an empty tuple"))
    S = eltype(children[1])
    isconcretetype(S) || throw(ArgumentError("property type must be concrete"))
    all(eltype(property) === S for property in children) ||
        throw(ArgumentError("all property children must have the same concrete type"))
    properties = Tuple(children)
    GeometryTupleProperties{S, typeof(properties)}(properties)
end

function get_value(property, ::Tuple)
    property
end

function get_value(properties::GeometryProperties, index::Tuple)
    properties isa GeometryLeafProperties && return properties.value
    isempty(index) && throw(ArgumentError("property index is empty for a nested property"))
    child_index = index[1]
    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    get_value(properties.properties[child_index], index[2:end])
end

function _get_value_many(properties::GeometryProperties, indices::Tuple)
    isempty(indices) && return 0
    properties isa GeometryLeafProperties && return properties.value

    grouped_indices = Dict{Int, Vector{Tuple}}()
    for index in indices
        isempty(index) && throw(ArgumentError("property indices must not be empty"))
        child_index = index[1]
        child_index isa Int || throw(ArgumentError("property child indices must be integers"))
        push!(get!(grouped_indices, child_index, Tuple[]), index[2:end])
    end

    sum(_get_value_many(
        properties.properties[child_index],
        Tuple(child_indices),
    ) for (child_index, child_indices) in grouped_indices)
end

function get_value(properties::GeometryProperties, indices::SVector)
    _get_value_many(properties, Tuple(indices))
end


function get_value(property, indices::AbstractVector{<:Tuple})
    isempty(indices) && return 0
    property
end

function get_value(properties::GeometryProperties, indices::AbstractVector{<:Tuple})
    isempty(indices) && return 0
    properties isa GeometryLeafProperties && return properties.value

    total = nothing
    child_index = nothing
    child_indices = Tuple[]
    for index in indices
        isempty(index) && throw(ArgumentError("inside indices must not be empty"))
        next_child_index, remaining = index[1], index[2:end]
        if child_index !== nothing && next_child_index < child_index
            throw(ArgumentError("inside indices must be sorted"))
        end
        if child_index !== nothing && next_child_index != child_index
            1 <= child_index <= length(properties.properties) ||
                throw(BoundsError(properties.properties, child_index))
            value = get_value(properties.properties[child_index], child_indices)
            total = total === nothing ? value : total + value
            empty!(child_indices)
        end
        child_index = next_child_index
        push!(child_indices, remaining)
    end

    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    value = get_value(properties.properties[child_index], child_indices)
    total === nothing ? value : total + value
end

function _all_property_values(properties::Tuple)
    all(property -> property isa GeometryLeafProperties, properties) && return (
        length(properties) == 1 ? properties[1].value : Tuple(property.value for property in properties)
        for _ in 1:1
    )

    if all(property -> property isa GeometryVectorProperties || property isa GeometryLeafProperties, properties)
        vector_properties = filter(property -> property isa GeometryVectorProperties, properties)
        length_properties = length(first(vector_properties).properties)
        all(length(property.properties) == length_properties for property in vector_properties) ||
            throw(ArgumentError("property vector structures must have matching lengths"))
        return Iterators.flatten(
            (
                _all_property_values(Tuple(
                    property isa GeometryLeafProperties ? property : property.properties[index]
                    for property in properties
                ))
                for index in 1:length_properties
            ),
        )
    end

    if all(property -> property isa GeometryTupleProperties || property isa GeometryLeafProperties, properties)
        tuple_properties = filter(property -> property isa GeometryTupleProperties, properties)
        length_properties = length(first(tuple_properties).properties)
        all(length(property.properties) == length_properties for property in tuple_properties) ||
            throw(ArgumentError("property tuple structures must have matching lengths"))
        return Iterators.flatten(
            (
                _all_property_values(Tuple(
                    property isa GeometryLeafProperties ? property : property.properties[index]
                    for property in properties
                ))
                for index in 1:length_properties
            ),
        )
    end

    throw(ArgumentError("property structures must match, apart from scalar leaves"))
end

"""
    all_property_values(properties...)

Iterate over all scalar values in one or more nested geometry properties. Scalar
leaves broadcast over nested properties, while non-leaf structures must match.
"""
function all_property_values(properties::GeometryProperties...)
    isempty(properties) && throw(ArgumentError("at least one property is required"))
    _all_property_values(properties)
end

end
