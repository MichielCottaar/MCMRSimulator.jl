module TimeSteps
import ..Constants: gyromagnetic_ratio
import ..Geometries.Internal: size_scale, max_timestep_sticking

"""
    TimeStep(simulation; timestep=Inf, turtoisity_precision=3e-2 * <precision>, gradient_precision=1e-4 * <precision>, precision=1.)

Creates an object controlling the timestep of the MCMR simulation.

At any time the timestep is guaranteed to be shorter than:
1. user-provided value of `timestep` (in ms).
2. `FullTimeStep.turtoisity_precision` * `size_scale(geometry)`^2 / D, where [`size_scale`](@ref) is the average size of the obstructions and `D` is the [`diffusivity`](@ref).
3. timestep that would allow permeability or magnetisation transfer probability to be close to 1.
4. (`FullTimeStep.gradient_precision` /( D * \\gamma^2 * G^2))^(1//3), where \\gamma is the [`gyromagnetic_ratio`](@ref) and `G` is the current `gradient_strength`.
"""
mutable struct TimeStep
    max_timestep :: Float64
    scaled_gradient_precision :: Float64
end


function TimeStep(; diffusivity, geometry, timestep=Inf, turtoisity_precision=nothing, gradient_precision=nothing, precision=1.)
    if iszero(diffusivity)
        return TimeStep(timestep, Inf)
    end
    return TimeStep(
        min(
            timestep,
            (isnothing(turtoisity_precision) ? (3e-2 * precision) : turtoisity_precision) * size_scale(geometry)^2 / diffusivity,
            max_timestep_sticking(geometry, diffusivity)
        ),
        (isnothing(gradient_precision) ? (1e-4 * precision) : gradient_precision) / diffusivity
    )
end

(ts::TimeStep)(gradient_strength) = min(
    ts.max_timestep,
    (ts.scaled_gradient_precision / gradient_strength^2)^(1//3)
)


end