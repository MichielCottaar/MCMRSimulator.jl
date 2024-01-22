module SequenceBuilders
import JuMP: Model
import Ipopt
import ..BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, match_blocks!
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

function SequenceBuilder(model::Model, blocks::AbstractVector) 
    SequenceBuilder(model, [to_block(model, b) for b in blocks])
end


function SequenceBuilder(blocks::AbstractVector)
    model = Model(Ipopt.Optimizer)
    SequenceBuilder(model, blocks)
end


end
