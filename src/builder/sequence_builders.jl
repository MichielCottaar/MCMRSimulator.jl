module SequenceBuilders
import JuMP: Model, owner_model, index, VariableRef, @constraint, @variable
import Ipopt
import ..BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, match_blocks!, duration
import ..Wait: WaitBlock

"""
    SequenceBuilder(blocks...)

Defines a sequence as a series of [`BuildingBlock`](@ref) objects.

After defining the blocks, the user can add one or more constraints and an objective function to the properties of the [`BuildingBlock`](@ref) objects.
A sequence matching these constraints will be produced by calling [`solve`](@ref)(builder) or [`Sequence`](@ref)(builder).
"""
struct SequenceBuilder
    model :: Model
    blocks :: Vector{<:BuildingBlock}
    TR :: VariableRef
end

function to_block(model::Model, placeholder::BuildingBlockPlaceholder{T}) where {T}
    block = T(model, placeholder.args...; placeholder.kwargs...)
    if isassigned(placeholder.concrete)
        match_blocks!(placeholder.concrete[], block)
    else
        placeholder.concrete[] = block
    end
    return block
end

to_block(model::Model, time::Union{Number, Symbol, Nothing, Val{:min}, Val{:max}}) = WaitBlock(model, time)

function SequenceBuilder(model::Model, blocks...) 
    builder = SequenceBuilder(model, BuildingBlock[], @variable(model))
    for b in blocks
        add!(builder.blocks, to_block(builder, b))
    end
    @constraint model TR(builder) >= duration(builder)
    return builder
end


function SequenceBuilder(blocks...)
    model = Model(Ipopt.Optimizer)
    SequenceBuilder(model, blocks)
end


function Base.show(io::IO, builder::SequenceBuilder)
    print(io, "SequenceBuilder($(builder.blocks)) being solved by ")
    show(io, builder.model)
end


# Interactions between individual BuildingBlock and parent
Base.length(sb::SequenceBuilder) = length(sb.blocks)
builder(bb::BuildingBlock) = bb.builder
owner_model(bb::BuildingBlock) = owner_model(builder(bb))
owner_model(sb::SequenceBuilder) = sb.model

"""
    TR(sequence::SequenceBuilder)

Return the repetition time of the sequence (in ms).
"""
TR(sb::SequenceBuilder) = sb.TR

"""
    index(bb::BuildingBlock)

Returns the index of the [`BuildingBlock`](@ref) in the parent [`SequenceBuilder`](@ref).
"""
index(sb::SequenceBuilder, bb::BuildingBlock) = findfirst(isequal(bb), sb.blocks)
index(bb::BuildingBlock) = index(builder(bb), bb)


"""
    duration(bb1::BuildingBlock, bb2::BuildingBlock)
    duration(sb::SequenceBuilder, bb1::Integer, bb2::Integer)

The duration of the sequence from the start of [`BuildingBlock`](@ref) `bb1` till the end of `bb2`.
Returns an error if they are not part of the same sequence.
If `bb1` plays after `bb2` then this will calculate the time between `bb2` and the `bb1` in the next repetition time
(i.e., `TR(sb) + end_time(bb2) - start_time(bb1)`)
"""
function duration(bb1::BuildingBlock, bb2::BuildingBlock)
    if builder(bb1) != builder(bb2)
        error("Cannot compute duration between blocks that are part of different sequences!")
    end
    sb = builder(bb1)
    return duration(sb, index(sb, bb1), index(sb, bb2))
end

function duration(sb::SequenceBuilder, index1::Integer, index2::Integer)
    if index2 == index1
        return duration(sb[index1])
    elseif index2 == index1 - 1
        return TR(sb)
    elseif index2 > index1
        return [duration(i) for i in index1:index2]
    else
        return [TR(sb) - duration(sb, index2+1, index1-1)]
    end
end

duration(sb::SequenceBuilder) = duration(sb, 1, length(sb))


"""
    start_time(building_block::BuildingBlock)
    start_time(seq::SequenceBuilder, building_block::Integer)

Return the start time of the given [`BuildingBlock`](@ref) within the [`SequenceBuilder`](@ref).
You can pass on the actual [`BuildingBlock`](@ref) object or the [`SequenceBuilder`](@ref) together with the index of the `building_block`
"""
start_time(bb::BuildingBlock) = start_time(builder(bb), bb)

function start_time(sb::SequenceBuilder, bb::BuildingBlock)
    if builder(bb) != sb
        error("No start time of $bb within $sb, because this BuildingBlock is not part of that sequence.")
    end
    return start_time(sb, index(bb))
end

start_time(sb::SequenceBuilder, index::Integer) = index == 1 ? 0. : duration(sb, 1, index - 1)

"""
    end_time(building_block::BuildingBlock)
    end_time(seq::SequenceBuilder, building_block::Integer)

Return the end time of the given [`BuildingBlock`](@ref) within the [`SequenceBuilder`](@ref).
You can pass on the actual [`BuildingBlock`](@ref) object or the [`SequenceBuilder`](@ref) together with the index of the `building_block`
"""
end_time(bb::BuildingBlock) = end_time(builder(bb), bb)

function end_time(sb::SequenceBuilder, bb::BuildingBlock)
    if builder(bb) != sb
        error("No end time of $bb within $sb, because this BuildingBlock is not part of that sequence.")
    end
    return end_time(sb, index(bb))
end

end_time(sb::SequenceBuilder, index::Integer) = duration(sb, 1, index)


end
