"""
Defines a set of different options for MRI gradients.
"""
module PulsedGradients

import JuMP: @constraint, @variable, Model, VariableRef
import StaticArrays: SVector
import ...BuildingBlocks: BuildingBlock, duration, helper_functions, set_simple_constraints!, scanner_constraints!
import ....Scanners: Scanner


"""
    PulsedGradient(; orientation=:bvec, variables...)

Defines a trapezoidal pulsed gradient

This is a [`BuildingBlock`](@ref) for a [`SequenceBuilder`](@ref).

## Parameters
- `orientation` sets the gradient orienation. Can be set to a vector for a fixed orientation. Alternatively, can be set to :bvec (default) to rotate with the user-provided `bvecs` or to :neg_bvec to always be the reverse of the `bvecs`.

## Variables
Variables can be set during construction or afterwards as an attribute.
If not set, they will be determined during the sequence optimisation.
### Timing variables
- [`rise_time`](@ref): Time of the gradient to reach from 0 to maximum in ms. If explicitly set to 0, the scanner slew rate will be ignored.
- [`flat_time`](@ref): Time that the gradient stays at maximum strength in ms.
- [`δ`](@ref): effective pulse duration (`rise_time` + `flat_time`) in ms.
- [`duration`](@ref): total pulse duration (2 * `rise_time` + `flat_time`) in ms.
### Gradient variables
- [`gradient_strength`](@ref): Maximum gradient strength achieved during the pulse in kHz/um
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um (given by `δ` * `gradient_strength`).

The [`bvalue`](@ref) can be constrained for multiple gradient pulses.
"""
mutable struct PulsedGradient <: BuildingBlock
    orientation :: Any
    slew_rate :: VariableRef
    rise_time :: VariableRef
    flat_time :: VariableRef
end

function PulsedGradient(model::Model; orientation=:bvec, kwargs...)
    res = PulsedGradient(
        orientation,
        @variable(model),
        @variable(model),
        @variable(model)
    )

    set_simple_constraints!(model, res, kwargs)
    @constraint model flat_time(res) >= 0
    @constraint model rise_time(res) >= 0
    @constraint model slew_rate(res) >= 0
    return res
end

"""
    rise_time(pulsed_gradient)

The time from 0 till the maximum gradient strength in ms.
"""
rise_time(pg::PulsedGradient) = pg.rise_time

"""The time spent at the maximum gradient strength in ms."""
flat_time(pg::PulsedGradient) = pg.flat_time

"""Maximum gradient strength in kHz/um."""
gradient_strength(g::PulsedGradient) = rise_time(g) * slew_rate(g)

"""Maximum rate of increase (and decrease) of the gradient strength in kHz/um/ms."""
slew_rate(g::PulsedGradient) = g.slew_rate

"""Pulse gradient duration (`rise_time + `flat_time).  This is the effective duration of the gradient. The real duration is longer (and given by [`duration`](@ref))."""
δ(g::PulsedGradient) = rise_time(g) + flat_time(g)

duration(g::PulsedGradient) = 2 * rise_time(g) + flat_time(g)

"""q-value at the end of the gradient (rad/um)."""
qval(g::PulsedGradient) = gradient_strength(g) * δ(g)

helper_functions(::Type{PulsedGradient}) = [qval, δ, gradient_strength, duration, rise_time, flat_time, slew_rate]

function scanner_constraints!(model::Model, g::PulsedGradient, scanner::Scanner)
    @constraint model gradient_strength(g) <= scanner.gradient
    @constraint model slew_rate(g) <= scanner.slew_rate
end

end