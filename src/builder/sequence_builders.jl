module SequenceBuilders
import JuMP: Model
import Ipopt
import ..BuildingBlocks: BuildingBlock
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


function to_block(block::Expr)
    if block.head == :call
        return Expr(:call, esc(block.args[1]), :model, esc.(block.args[2:end])...)
    else
        error("$(block) is not a valid part of a sequence.")
    end
end

function to_block(block::Union{QuoteNode, Number, Symbol})
    return :(WaitBlock(model, $(esc(block)))) 
end

macro builder(blocks...)
    model_blocks = to_block.(blocks)
    create_builder = Expr(:call, :SequenceBuilder, :model, Expr(:vect, model_blocks...))
    quote
        model = Model(Ipopt.Optimizer)
        $(create_builder)
    end
end

function Base.show(io::IO, builder::SequenceBuilder)
    print(io, "SequenceBuilder(")
    for block in builder.blocks
        print(io, block, ", ")
    end
    print(io, ")")
end

end
