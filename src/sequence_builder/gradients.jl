"""
Defines a set of different options for MRI gradients.
"""
module Gradients
import JuMP: VariableRef, @constraint, @variable, Model
import StaticArrays: SVector
import ..BuildingBlocks: BuildingBlock, duration, helper_functions, show_helper_functions, apply_scanner!, set_kwargs_constraints!, add_parameters_to_docs
import ...Scanners: Scanner


struct PulsedGradient <: BuildingBlock
    slew_rate :: VariableRef
    flat_time :: VariableRef
    rise_time :: VariableRef
    orientation :: Union{Nothing, SVector{3, Float64}}
end

function PulsedGradient(model::Model; orientation=nothing, kwargs...)
    if !isnothing(orientation)
        orientation = SVector{3, Float64}(orientation)
        if all(iszero.(orientation))
            error("Gradient orienation is invalid. At least one of the values should be non-zero.")
        end
    end
    res = PulsedGradient(
        @variable(model),
        @variable(model),
        @variable(model),
        orientation
    )
    @constraint model rise_time(res) >= 0.
    @constraint model flat_time(res) >= 0.
    @constraint model slew_rate(res) >= 0.

    set_kwargs_constraints!(model, res; kwargs...)
    return res
end

PulsedGradient(; kwargs...) = (PulsedGradient, kwargs)


"""
    rise_time(pulsed_gradient)

The time from 0 till the maximum gradient strength in ms.
"""
rise_time(g::PulsedGradient) = g.rise_time

"""
    flat_time(pulsed_gradient)

The time spent at the maximum gradient strength in ms.
"""
flat_time(g::PulsedGradient) = g.flat_time

"""
    gradient_strength(gradient)

Maximum gradient strength in kHz/um.
"""
gradient_strength(g::PulsedGradient) = rise_time(g) * slew_rate(g)

"""
    slew_rate(gradient)

Maximum rate of increase (and decrease) of the gradient strength in kHz/um/ms.
"""
slew_rate(g::PulsedGradient) =  g.slew_rate

"""
    δ(pulsed_gradient)

Pulse gradient duration ([`rise_time`](@ref) + [`flat_time`](@ref)).
This is the effective duration of the gradient.
The real duration is longer (and given by [`duration`](@ref)).
"""
δ(g::PulsedGradient) = rise_time(g) + flat_time(g)

duration(g::PulsedGradient) = 2 * rise_time(g) + flat_time(g)

"""
    qval(gradient)

q-value at the end of the gradient.
"""
qval(g::PulsedGradient) = gradient_strength(g) * δ(g)

helper_functions(::Type{PulsedGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time, slew_rate]
show_helper_functions(::Type{PulsedGradient}) = [(δ, "ms"), (gradient_strength, "kHz/um"), (qval, "kHz"), (rise_time, "ms")]
@doc add_parameters_to_docs(PulsedGradient, """
    PulsedGradient(parameters...)

A [`PulsedGradient`](@ref) represents a single gradient pulse (i.e., a trapezoidal gradient profile).

## Parameters
- `orientation`: If a vector is provided, the gradient will always point in the direction of that vector. Otherwise, it will rotate with the user-supplied `bvecs` (and point in the x-direction if no `bvecs` are supplied).
""") PulsedGradient

function apply_scanner!(model::Model, scanner::Scanner, g::PulsedGradient)
    if isnothing(g.orientation)
        directional_excess = 1
    else
        directional_excess = norm(g.orientation) / maximum(abs.(g.orientation))
    end
    @constraint model gradient_strength(g) <= scanner.gradient * directional_excess
    @constraint model slew_rate(g) <= scanner.slew_rate * directional_excess
end



end