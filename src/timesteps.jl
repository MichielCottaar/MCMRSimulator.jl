module TimeSteps
import ..Constants: gyromagnetic_ratio
import ..Geometries: Internal

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


function TimeStep(; diffusivity, geometry, timestep=Inf, size_scale=nothing, turtoisity_precision=nothing, gradient_precision=nothing, precision=1., verbose=true)
    if iszero(diffusivity)
        return TimeStep(timestep, Inf)
    end
    use_size_scale = isnothing(size_scale) ? Internal.size_scale(geometry) : size_scale
    options = (
            timestep,
            (isnothing(turtoisity_precision) ? (3e-2 * precision) : turtoisity_precision) * use_size_scale^2 / diffusivity,
            Internal.max_timestep_sticking(geometry, diffusivity)
        )
    idx = argmin(options)
    if verbose
        @info "# Timestep determination"
        if idx == 1
            @info "Maximum timestep set by user to $(options[1]) ms."
        elseif idx == 2
            @info "Maximum timestep set by turtoisity constraint based on size of geometry to $(options[2]) ms."
            if isnothing(size_scale)
                @info "Size scale of smallest object in the simulation was automatically determined to be $(use_size_scale) um."
                @info "If this value is too small, you can set the size scale explicitly by passing on `size_scale=<new_value>` to the `Simulator` constructor."
            else
                @info "Size scale was set to $(use_size_scale) um by the user."
            end
        elseif idx == 3
            @info "Maximum timestep set by requirement to get sufficient transitions from free to bound spins to $(options[3]) ms."
        end
        @info "The actual timestep might be further reduced based on the MR sequence(s)."
    end
    return TimeStep(minimum(options), (isnothing(gradient_precision) ? (1e-4 * precision) : gradient_precision) / diffusivity)
end

(ts::TimeStep)(gradient_strength) = min(
    ts.max_timestep,
    (ts.scaled_gradient_precision / gradient_strength^2)^(1//3)
)


end