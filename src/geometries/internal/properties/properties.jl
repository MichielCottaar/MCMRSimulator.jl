"""MRI-property storage and lookup for physical geometries."""
module Properties

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

"""Return the property value and the path consumed before reaching its scalar leaf."""
function get_value(property, ::Tuple)
    property, ()
end

get_value(properties::GeometryLeafProperties, ::Tuple) = (properties.value, ())

function get_value(properties::GeometryProperties, index::Tuple)
    isempty(index) && throw(ArgumentError("property index is empty for a nested property"))
    child_index = index[1]
    value, consumed = get_value(properties.properties[child_index], Base.tail(index))
    value, (child_index, consumed...)
end

"""Sum each distinct property-tree leaf once, assuming that the `indices` are ordered."""
function sum_values(properties::GeometryProperties, indices)
    total = zero(eltype(properties))
    previous = nothing
    for index in indices
        value, consumed = get_value(properties, index)
        if consumed == previous
            continue
        end
        previous = consumed
        total += value
    end
    total
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
