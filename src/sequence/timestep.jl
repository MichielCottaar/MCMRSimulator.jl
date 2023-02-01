struct TimeController
    timestep :: Float
    gradient_precision :: Float
    sample_displacement :: Int
    sample_off_resonance :: Int
    rf_rotation :: Float
end

TimeController(timestep :: Number) = TimeController(Float(timestep), 0., 0, 0, 0.)
function TimeController(geometry; gradient_precision=0.01, sample_displacement=5, sample_off_resonance=10, rf_rotation=1.)
    if length(geometry) == 0
        sample_displacement = 1
    end
    if !produces_off_resonance(geometry)
        sample_off_resonance = 1
    end
    TimeController(zero(Float), Float(gradient_precision), Int(sample_displacement), Int(sample_off_resonance), Float(rf_rotation))
end

"""
    propose_times(simulation::Simulation, t_start, t_end)
    propose_times(time_controller::TimeController, t_start, t_end, sequences, diffusivity)

Computes the timepoints at which the simulation will be evaluated when running from `t_start` to `t_end`.

The following timepoints will always be included
- Any multiple of the TR
- Any RF pulse
- Any change in the gradient strength

Additional timepoints will be added to ensure that:
- at any time the timestep is larger than `time_controller.gradient_precision` times ``D^{-1/3} (G \\gamma)^{-2/3}``.
- between `t_start` and `t_end` there are at least `time_controller.sample_displacement` timepoints.
- between any control point of the sequence gradient there are at least `time_controller.sample_displacement` timepoints.
- between an RF pulse and a subsequent RF pulse or t_end there are at least `time_controller.sample_off_resonance` timepoints.
- During an RF block the rotation around the maximum magnetic field will be at most `time_controller.rf_rotation` radians. For a specific [`RFBlock`](@ref) these timepoints can be found using `propose_times(rf_block, rf_rotation, B0)`.
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
        Ntimepoints = Int(div(t1 - t0, max_timestep(sequences, time_controller, t0, t1, t_end, diffusivity), RoundUp))
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
        control_points(sequence.gradient)...,
        start_time.(sequence.pulses)...,
        end_time.(sequence.pulses)...,
    ]


function max_timestep(sequences::AbstractVector{<:Sequence}, time_controller::TimeController, t_start::Number, t_end::Number, next_readout::Number, diffusivity::Number)
    if !iszero(time_controller.timestep)
        return time_controller.timestep
    elseif length(sequences) == 0
        return max_timestep_displacement_readouts(time_controller.sample_displacement, t_start, t_end)
    else
        return min(
            max_timestep_gradient(sequences, time_controller.gradient_precision, t_start, t_end, diffusivity),
            max_timestep_displacement(sequences, time_controller.sample_displacement, t_start, t_end),
            max_timestep_displacement_readouts(time_controller.sample_displacement, t_start, t_end),
            max_timestep_off_resonance(sequences, time_controller.sample_off_resonance, t_start, t_end, next_readout),
            max_timestep_pulse(sequences, time_controller.rf_rotation, t_start, t_end),
        )
    end
end

for string in ("gradient", "displacement", "off_resonance", "pulse")
    symb = Symbol("max_timestep_" * string)
    @eval $symb(sequences::AbstractVector{<:Sequence}, args...) = minimum([$symb(s, args...) for s in sequences])
end


function max_timestep_gradient(sequence::Sequence, gradient_precision::Number, t_start::Number, t_end::Number, diffusivity::Number)
    if iszero(gradient_precision)
        return Inf
    end
    grad = norm(gradient(sequence, t_start, t_end))  # frequency gradient in kHz/um
    return gradient_precision * (diffusivity * grad * grad) ^ (-1/3)
end


max_timestep_displacement_readouts(sample_displacement::Number, t_start::Number, t_end::Number) = iszero(sample_displacement) ? Inf : (t_end - t_start) / sample_displacement

function max_timestep_displacement(sequence::Sequence, sample_displacement::Number, t_start::Number, t_end::Number)
    if iszero(sample_displacement)
        return Inf
    end
    mean_time = (t_start + t_end) / 2
    cp = control_points(sequence.gradient)
    index_control = searchsortedfirst(cp, mod(mean_time, sequence.TR))
    if index_control <= 1 || (index_control >= length(cp) + 1)
        return Inf
    end
    return (
        cp[index_control] -
        cp[index_control - 1]
    ) / sample_displacement
end
    
function max_timestep_off_resonance(sequence::Sequence, sample_off_resonance::Number, t_start::Number, t_end::Number, next_readout::Number)
    if iszero(sample_off_resonance)
        return Inf
    end
    mean_time = (t_start + t_end) / 2
    if !isnothing(current_pulse(sequence, mean_time))
        # ignore this if we are currently in an RF pulse
        return Inf
    end
    pi = previous_instant(sequence, mean_time)
    pp = previous_pulse(sequence, mean_time)
    if isnothing(pi) && isnothing(pp)
        # We are before the first RF pulse
        return Inf
    end
    start_time_default(n) = isnothing(n) ? Inf : start_time(n)
    next_timepoint = min(
        next_readout, 
        start_time_default(next_instant(sequence, mean_time)),
        start_time_default(next_pulse(sequence, mean_time)),
    )
    end_time_default(n) = isnothing(n) ? 0. : end_time(n)
    prev_timepoint = max(
        end_time_default(pi),
        end_time_default(pp),
    )
    return (next_timepoint - prev_timepoint) / sample_off_resonance
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