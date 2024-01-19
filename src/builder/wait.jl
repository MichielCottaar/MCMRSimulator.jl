module Wait
import JuMP: Model, @constraint, @variable, VariableRef
import ..BuildingBlocks: BuildingBlock, duration, helper_functions, apply_simple_constraint!
import ...Scanners: Scanner

"""
    WaitBlock(duration)

An empty [`BuildingBlock`](@ref) of given `duration` (in ms).

Duration can be set to one of:
- numeric value to fix it
- `:min` to minimise its value given any external constraints
- `:max` to maximise its value given any external constraints
- `nothing` to make it fully determined by external constraints and objectives
"""
struct WaitBlock <: BuildingBlock
    duration :: VariableRef
end

function WaitBlock(model::Model, duration_constraint=nothing)
    res = WaitBlock(@variable(model))
    @constraint model duration(res) >= 0
    if !isnothing(duration_constraint)
        apply_simple_constraint!(model, duration(res), duration_constraint)
    end
    return res
end

helper_functions(::Type{WaitBlock}) = [duration]

duration(wb::WaitBlock) = wb.duration

scanner_constraints!(::Model, ::WaitBlock, ::Scanner) = nothing

end