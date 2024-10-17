module TimeSteps
import ..Constants: gyromagnetic_ratio
import ..Geometries: Internal

"""
    TimeStep(simulation; turtoisity=3e-2, gradient=1e-4, permeability=0.5,)

Creates an object controlling the timestep of the MCMR simulation.

At any time the timestep is guaranteed to be shorter than:
1. `FullTimeStep.turtoisity_precision` * `size_scale(geometry)`^2 / D, where [`size_scale`](@ref) is the average size of the obstructions and `D` is the [`diffusivity`](@ref).
2. timestep greater than `permeability` times 1 / (maximum permeability parameter)^2
3. timestep that would allow magnetisation transfer probability to be greater than 25%.
4. the minimum dwell time of the bound pool.
5. (`FullTimeStep.gradient_precision` /( D * \\gamma^2 * G^2))^(1//3), where \\gamma is the [`gyromagnetic_ratio`](@ref) and `G` is the current `gradient_strength`.
"""
mutable struct TimeStep
    max_timestep :: Float64
    scaled_gradient_precision :: Float64
end


function TimeStep(; 
    diffusivity, geometry, size_scale=nothing, 
    turtoisity=3e-2, gradient=1e-4, verbose=true,
    max_permeability_probability=0.5,
    )
    if iszero(diffusivity)
        return TimeStep(Inf, Inf)
    end
    use_size_scale = isnothing(size_scale) ? Internal.size_scale(geometry) : size_scale
    options = (
            timestep,
            turtoisity * use_size_scale^2 / diffusivity,
            Internal.max_timestep_sticking(geometry, diffusivity),
            max_timestep_permeability(geometry, max_permeability_probability),
        )
    idx = argmin(options)
    if verbose
        lines = ["# Timestep determination"]
        if idx == 1
            push!(lines, "Maximum timestep set by user to $(options[1]) ms.")
        elseif idx == 2
            push!(lines, "Maximum timestep set by turtoisity constraint based on size of geometry to $(options[2]) ms.")
            if isnothing(size_scale)
                push!(lines, "Size scale of smallest object in the simulation was automatically determined to be $(use_size_scale) um.")
                push!(lines, "If this value is too small, you can set the size scale explicitly by passing on `size_scale=<new_value>` to the `Simulator` constructor.")
            else
                push!(lines, "Size scale was set to $(use_size_scale) um by the user.")
            end
        elseif idx == 3
            push!(lines, "Maximum timestep set by requirement to get sufficient transitions from free to bound spins to $(options[3]) ms.")
        elseif idx == 4
            push!(lines, "Maximum timestep set by requirement to get sufficient rate of spins through the permeable surfaces to $(options[4]) ms.")
        end
        push!(lines, "The actual timestep might be further reduced based on the MR sequence(s).")
        @info join(lines, '\n')
    end
    return TimeStep(minimum(options), gradient / diffusivity)
end

(ts::TimeStep)(gradient_strength) = min(
    ts.max_timestep,
    (ts.scaled_gradient_precision / gradient_strength^2)^(1//3)
)

"""
    max_timestep_permeability(geometry, scaling)

Compute the maximum timestep necessary to keep the probability of a spin passing through sufficiently low.
"""
max_timestep_permeability(geometry, scaling) = scaling * (1. / Internal.max_permeability_non_inf(geometry))^2


end