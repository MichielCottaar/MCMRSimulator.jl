"""MRI-property storage and lookup for physical geometries."""
module Properties

import StaticArrays: SVector
import ..Indices: ObstructionIndex

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

function get_value(property, ::ObstructionIndex)
    property
end

function get_value(properties::GeometryProperties, index::ObstructionIndex)
    properties isa GeometryLeafProperties && return properties.value
    child_index = index.indices[1]
    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    remaining_length = length(index.indices) - 1
    remaining = ObstructionIndex{remaining_length}(
        SVector{remaining_length, Int}(index.indices[2:end]),
    )
    get_value(properties.properties[child_index], remaining)
end


function get_value(property, indices::AbstractVector{<:ObstructionIndex})
    isempty(indices) && return 0
    property
end

function get_value(properties::GeometryProperties, indices::AbstractVector{<:ObstructionIndex})
    isempty(indices) && return 0
    properties isa GeometryLeafProperties && return properties.value

    total = nothing
    child_index = nothing
    child_indices = ObstructionIndex[]
    for index in indices
        isempty(index.indices) && throw(ArgumentError("inside indices must not be empty"))
        next_child_index = index.indices[1]
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
        remaining_length = length(index.indices) - 1
        push!(child_indices, ObstructionIndex{remaining_length}(
            SVector{remaining_length, Int}(index.indices[2:end]),
        ))
    end

    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    value = get_value(properties.properties[child_index], child_indices)
    total === nothing ? value : total + value
end

end
