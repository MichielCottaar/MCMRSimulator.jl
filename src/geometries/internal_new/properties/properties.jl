"""MRI-property storage and lookup for physical geometries."""
module Properties

import StaticArrays: SVector
import ..Indices: ObstructionIndex

struct GeometryVectorProperties{P}
    properties::Vector{P}
end

struct GeometryTupleProperties{P<:Tuple}
    properties::P
end

function get_value(property, ::ObstructionIndex)
    property
end

function get_value(properties::GeometryVectorProperties, index::ObstructionIndex)
    child_index = index.indices[1]
    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    remaining_length = length(index.indices) - 1
    remaining = ObstructionIndex{remaining_length}(
        SVector{remaining_length, Int}(index.indices[2:end]),
    )
    get_value(properties.properties[child_index], remaining)
end

function get_value(properties::GeometryTupleProperties, index::ObstructionIndex)
    child_index = index.indices[1]
    1 <= child_index <= length(properties.properties) ||
        throw(BoundsError(properties.properties, child_index))
    remaining_length = length(index.indices) - 1
    remaining = ObstructionIndex{remaining_length}(
        SVector{remaining_length, Int}(index.indices[2:end]),
    )
    get_value(properties.properties[child_index], remaining)
end

get_value(property, indices::Vector{<:ObstructionIndex}) = [
    get_value(property, index) for index in indices
]

end
