module InstantPulses
import JuMP: @constraint, @variable, VariableRef
import ...BuildingBlocks: BuildingBlock, properties, BuildingBlockPlaceholder, set_simple_constraints!, duration
import ...SequenceBuilders: SequenceBuilder, owner_model

struct InstantRFPulseBlock <: BuildingBlock
    builder :: SequenceBuilder
    flip_angle :: VariableRef
    phase :: VariableRef
end

InstantRFPulseBlock(; kwargs...) = BuildingBlockPlaceholder{InstantRFPulseBlock}(; kwargs...)
function InstantRFPulseBlock(builder::SequenceBuilder; kwargs...) 
    model = owner_model(builder)
    res = InstantRFPulseBlock(
        builder,
        @variable(model),
        @variable(model)
    )
    @constraint model flip_angle(res) >= 0
    set_simple_constraints!(model, res, kwargs)
    return res
end

flip_angle(instant::InstantRFPulseBlock) = instant.flip_angle
phase(instant::InstantRFPulseBlock) = instant.phase
duration(instant::InstantRFPulseBlock) = 0.
properties(::Type{<:InstantRFPulseBlock}) = [flip_angle, phase]

end