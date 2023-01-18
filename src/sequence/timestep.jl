struct TimeController
    timestep :: Float
    gradient_precision :: Float
    sample_displacement :: Int
    sample_off_resonance :: Int
end

TimeController(timestep :: Float) = TimeController(timestep, 0., 0, 0)
function TimeController(geometry; gradient_precision=0.01, sample_displacement=5, sample_off_resonance=10)
    if length(geometry) == 0
        sample_displacement = 1
    end
    if !produces_off_resonance(geometry)
        sample_off_resonance = 1
    end
    TimeController(0., gradient_precision, sample_displacement, sample_off_resonance)
end

"""
    get_times(simulation::Simulation, t_start, t_end)
    get_times(time_controller::TimeController, t_start, t_end, sequences, diffusivity)

Computes the timepoints at which the simulation will be evaluated when running from `tstart` to `tend`.

The following timepoints will always be included
- Any multiple of the TR
- Any RF pulse
- Any change in the gradient strength

Additional timepoints will be added to ensure that:
- at any time the timestep is larger than `time_controller.gradient_precision` times ``D^{-1/3} (G \\gamma)^{-2/3}``.
- between `t_start` and `t_end` there are at least `time_controller.sample_displacement` timepoints.
- between any control point of the sequence gradient there are at least `time_controller.sample_displacement` timepoints.
- between an RF pulse and a subsequent RF pulse or t_end there are at least `time_controller.sample_off_resonance` timepoints.
"""
function get_times(time_controller::TimeController, t_start::Float, t_end::Float, sequences::AbstractVector{<:Sequence}, diffusivity::Float)
    timepoints = [t_start, t_end]
    for sequence in sequences
        index = 1
        while time(sequence, index) <= t_start
            index += 1
        end
        while time(sequence, index) < t_end
            push!(timepoints, time(sequence, index))
            index += 1
        end

        nTR = div(t_start, sequence.TR, RoundDown)
        all_control_points = sequence.gradient.times .+ (nTR * sequence.TR)
        control_points = all_control_points[all_control_points .> t_start]
        while (nTR + 1) * sequence.TR < t_end
            # This TR is fully included before `t_end`
            append!(timepoints, control_points)
            nTR += 1
            control_points = sequence.gradient.times .+ (nTR * sequence.TR)
        end
        append!(timepoints, control_points[control_points .< t_end])
    end
    timepoints = sort(unique(timepoints))

    if !iszero(time_controller.timestep)
        # Use the user-provided timestep
        for (t0, t1) in zip(copy(timepoints[1:end-1]), copy(timepoints[2:end]))
            Ntimepoints = div(t1 - t0, time_controller.timestep, RoundUp)
            if Ntimepoints > 1
                for new_timestep in range(t0, t1, length=Ntimesteps+1)[2:end-1]
                    push!(timepoints, new_timestep)
                end
            end
        end
        return sort(unique(timepoints))
    end

    for (t0, t1) in zip(copy(timepoints[1:end-1]), copy(timepoints[2:end]))
        mean_time = (t0 + t1) / 2

        # gradient precision
        max_timestep_gradient = Inf64
        if !iszero(time_controller.gradient_precision)
            for sequence in sequences
                grad = norm(get_gradient(sequence, t0, t1)) * gyromagnetic_ratio * Float(1e-3) * sequence.scanner.B0  # frequency gradient in rad / (ms um)
                mt = time_controller.gradient_precision * (diffusivity * grad * grad) ^ (-1/3)
                if mt < max_timestep_gradient
                    max_timestep_gradient = mt
                end
            end
        end

        # sample displacement between any two readouts and any two changes in gradient profile
        max_timestep_displacement = (t_end - t_start) / time_controller.sample_displacement
        for sequence in sequences
            control_points = sequence.gradient.times
            index_control = searchsortedfirst(control_points, mean_time)
            mt = (
                control_points[index_control] -
                control_points[index_control - 1]
            ) / time_controller.sample_frequency
            if mt < max_timestep_displacement
                max_timestep_displacement = mt
            end
        end

        # sample off-resonance field between RF pulse and subsequence RF pulse or readout
        max_timestep_offresonance = Inf64
        for sequence in sequences
            next_rf_pulse = next_pulse(sequence, mean_time)
            if next_rf_pulse <= 1
                continue
            end

            mt = (
                min(t_end, time(sequence, next_rf_pulse)) - 
                time(sequence, next_rf_pulse - 1)
            )
            if mt < max_timestep_offresonance
                max_timestep_offresonance = mt
            end
        end

        max_timestep = min(max_timestep_gradient, max_timestep_displacement, max_timestep_offresonance)

        if ~isinf(max_timestep)
            Ntimesteps = Int(div(t1 - t0, max_timestep, RoundUp))
            for new_timestep in range(t0, t1, length=Ntimesteps+1)[2:end-1]
                push!(timepoints, new_timestep)
            end
        end
    end
    sort!(timepoints)
    return timepoints
end
