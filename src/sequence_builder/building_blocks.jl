module BuildingBlocks
import ...Sequences: RFPulse, MRGradients, InstantComponent, Readout, Sequence, start_time, end_time, add_TR
import ...Scanners: Scanner

"""
    BuildingBlock(; components=vector of pulses/gradients, duration=minimum)

Creates a sequence building block by overlapping zero or more pulses/gradients.
"""
struct BuildingBlock
    components :: Vector
    duration :: Float64
    function BuildingBlock(; components=RFPulse[], duration=nothing)
        min_duration = length(components) == 0 ? 0 : maximum(end_time.(components))
        if isnothing(duration)
            duration = min_duration
        end
        @assert duration >= min_duration
        new(components, Float64(duration))
    end
end

const BuildingBlockLike = Union{BuildingBlock, RFPulse, MRGradients, Number, InstantComponent, Readout}

"""
    BuildingBlock(delay::Number)

Creates an empty sequence building block representing a delay of `delay` ms.
"""
BuildingBlock(number::Number) = BuildingBlock(duration=abs(number) < 1e-8 ? 0 : number)

"""
    BuildingBlock(pulse/gradient)

Creates an building block containing just the [`RFPulse`](@ref), [`InstantRFPulse`](@ref), [`MRGradients`](@ref), [`InstantGradient`](@ref), or [`Readout`].
The duration will be set to the length of the pulse/gradient.
"""
BuildingBlock(pulse::Union{RFPulse, MRGradients, InstantComponent, Readout}) = BuildingBlock(components=[pulse])
BuildingBlock(block::BuildingBlock) = block

"""
    BuildingBlock([building block like objects...])

Creates a sequence building block by concatenating all the [`BuildingBlock`](@ref)-like objects in the vector.
Each building block will play out sequenctially.
"""
function BuildingBlock(blocks::AbstractVector)
    all_components = []
    current_time = 0.
    for raw_block in blocks
        block = BuildingBlock(raw_block)
        append!(all_components, [add_TR(c, current_time) for c in block.components])
        current_time += block.duration
    end
    return BuildingBlock(components=all_components, duration=current_time)
end


function Sequence(block_like; TR=nothing, scanner=Scanner(B0=3.)) 
    block = BuildingBlock(block_like)
    Sequence(
        TR=isnothing(TR) ? block.duration : TR, 
        components=block.components,
        scanner=scanner
    )
end




"""
    duration(building_block)

Returns how long each sequence [`BuildingBlock`](@ref) lasts.
"""
duration(blocks::AbstractVector) = sum(duration.(blocks))
duration(block::BuildingBlockLike) = end_time(block) - start_time(block)
duration(block::BuildingBlock) = block.duration
duration(number::Number) = number
duration(r::Readout) = 0.


"""
    isempty_block(building_block)

Returns whether a sequence [`BuildingBlock`](@ref) is empty (i.e., does not contain RF pulses, readouts or MRI gradients).
"""
isempty_block(blocks::AbstractVector) = all(is_empty.(blocks))
isempty_block(number::Number) = true
isempty_block(block::BuildingBlock) = (length(block.pulses) == 0) && (length(block.instants) == 0) && (length(block.gradients) == 0)
isempty_block(other::BuildingBlockLike) = false
end