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

function get_value(properties::Union{GeometryVectorProperties, GeometryTupleProperties}, index::ObstructionIndex)
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

function get_value(
    properties::Union{GeometryVectorProperties, GeometryTupleProperties},
    indices::AbstractVector{<:ObstructionIndex},
)
    isempty(indices) && return 0

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
