"""
    TimeController(geometry, B0, diffusivity; max_stepsize=<see docs>, gradient_precision=1, rf_rotation=1)

Stores the settings controlling the maximum timestep during the simulation, namely:
- `max_timestep`: Generic maximum timestep that is considered throughout the simulation. By default this is set by the minimum of the following three options:
    1. ``\\Delta t_1 = l^2/(2 D)``, where `l` is the [`size_scale`](@ref) of the geometry, and ``D`` is the `diffusivity`.
    2. ``\\Delta t_2 = D^{-1/3} (360 * G)^{-2/3} / p``, where ``p`` is the `gradient_precision`, and ``G`` is the maximum off-resonance gradient due to the geometry (see [`off_resonance_gradient`](@ref)).
    3. The timestep required to keep the probability of the spins getting stuck on the surface at a collision below 1.
- `gradient_precision`: maximum error in the phase that a spin should incur while a gradient is active in any of the sequences (units: degrees). 
    If `max_stepsize` is not set by the user, it is also applied to any geometry-applied gradients (see point 2 above).
- `rf_rotation`: Maximum amount of rotation that an RF pulse should be able to apply in a single timestep (units: degrees).
    This constraint is only active if there is an RF pulse active in any of the sequences.

More details about how these settings are used can be found in [`propose_times`](ref).
"""
struct TimeController
    max_timestep :: Float
    gradient_precision :: Float
    rf_rotation :: Float
end

function TimeController(geometry::Geometry, B0::Number, diffusivity::Number, default_properties::GlobalProperties; max_timestep=nothing, gradient_precision=1., rf_rotation=1.)
    mt_stick = max_timestep_sticking(geometry, default_properties, diffusivity)
    if isnothing(max_timestep)
        max_timestep = min(
            size_scale(geometry),
            max_timestep_internal_gradient(geometry, gradient_precision, diffusivity, B0),
            mt_stick
        )
    elseif max_timestep > mt_stick
        throw(DomainError("Maximum timestep is set by the user to $(max_timestep)ms. However, it can not be longer than $(mt_stick)ms, because the probability of particles sticking to the surface would be larger than 1."))
    end
    TimeController(Float(max_timestep), Float(gradient_precision), Float(rf_rotation))
end

"""
    propose_times(simulation, t_start, t_end)
    propose_times(time_controller, t_start, t_end, sequences, diffusivity)

Computes the timepoints at which the simulation will be evaluated when running from `t_start` to `t_end`.

The following timepoints will always be included
- Any multiple of the TR
- Any RF pulse
- Any change in the gradient strength

Additional timepoints will be added to ensure that at any step the timestep is lower than the maximum timestep as described in [`TimeController`](@ref):
- the timestep needed to displace with the `simulation.time_controller.max_stepsize`.
- at any time the timestep is larger than ``D^{-1/3} (360 * G \\gamma)^{-2/3} / p``, where ``p`` is the `simulation.time_controller.gradient_precision` (given in degrees).
    The gradient ``G`` is the maximum of the user-applied gradient or the maximum off-resonance field gradient.
- During an [`RFPulse`](@ref) the rotation around the maximum magnetic field will be at most `simulation.time_controller.rf_rotation` degrees. 
"""
function propose_times(time_controller::TimeController, t_start::Number, t_end::Number, sequences::AbstractVector{<:Sequence}, diffusivity::Number)
    timepoints = Float[t_start, t_end]
    for sequence in sequences
        first_TR = all_control_points(sequence)

        min_nTR = div(t_start, sequence.TR, RoundDown)
        max_nTR = div(t_end, sequence.TR, RoundUp)
        for nTR in min_nTR:max_nTR
            TR_control_points = first_TR .+ (nTR * sequence.TR)
            control_points = TR_control_points[
                (TR_control_points .> t_start) .&
                (TR_control_points .< t_end)
            ]
            append!(timepoints, control_points)
        end
    end
    timepoints = sort(unique(timepoints))

    for (t0, t1) in zip(copy(timepoints[1:end-1]), copy(timepoints[2:end]))
        mt = max_timestep(sequences, time_controller, t0, t1, diffusivity)
        if isinf(mt)
            continue
        end
        Ntimepoints = Int(div(t1 - t0, mt, RoundUp))
        if Ntimepoints > 1
            for new_timepoint in range(t0, t1, length=Ntimepoints+1)[2:end-1]
                push!(timepoints, new_timepoint)
            end
        end
    end
    return sort(unique(timepoints))
end


all_control_points(sequence::Sequence) = [
        0, sequence.TR,
        [get_time(instant) for instant in sequence.instants]...,
        vcat(control_points.(sequence.gradients)...)...,
        start_time.(sequence.pulses)...,
        end_time.(sequence.pulses)...,
    ]


function max_timestep(sequences::AbstractVector{<:Sequence}, time_controller::TimeController, t_start::Number, t_end::Number, diffusivity::Number)
    if iszero(length(sequences))
        return time_controller.max_timestep
    else
        return min(
            time_controller.max_timestep,
            max_timestep_sequence_gradient(sequences, time_controller.gradient_precision, t_start, t_end, diffusivity),
            max_timestep_pulse(sequences, time_controller.rf_rotation, t_start, t_end),
        )
    end
end

for string in ("sequence_gradient", "pulse")
    symb = Symbol("max_timestep_" * string)
    @eval $symb(sequences::AbstractVector{<:Sequence}, args...) = minimum([$symb(s, args...) for s in sequences])
end

function max_timestep_sequence_gradient(sequence::Sequence, gradient_precision::Number, t_start::Number, t_end::Number, diffusivity::Number)
    if iszero(gradient_precision)
        return Inf
    end
    grad = norm(gradient(sequence, (t_start + t_end) / 2))
    return (360 * diffusivity * grad * grad) ^ (-1/3) / gradient_precision
end

function max_timestep_internal_gradient(geom::Geometry, gradient_precision::Number, diffusivity::Number, B0)
    if iszero(gradient_precision)
        return Inf
    end
    grad = off_resonance_gradient(geom, B0)
    return (360 * diffusivity * grad * grad) ^ (-1/3) / gradient_precision
end

function max_timestep_pulse(sequence::Sequence, rf_rotation::Number, t_start::Number, t_end::Number)
    if iszero(rf_rotation)
        return Inf
    end
    current = current_pulse(sequence, (t_start + t_end) / 2)
    if isnothing(current)
        return Inf
    end
    return rf_rotation / (current.max_amplitude * 360)
end