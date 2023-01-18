struct TimeController
    timestep :: Float
    gradient_precision :: Float
    sample_displacement :: Int
    sample_off_resonance :: Int
end

TimeController(timestep :: Float) = TimeController(timestep, 0., 0, 0)
function TimeController(simulation; gradient_prection=0.01, sample_displacement=5, sample_off_resonance=10)
    if length(simulation.micro.geometry) == 0
        sample_displacement = 1
    end
    if !produces_off_resonance(simulation)
        sample_off_resonance = 1
    end
    TimeController(0., gradient_prection, sample_displacement, sample_off_resonance)
end

"""
    get_times(simulation::Simulation, readout_times)
    get_times(time_controller::TimeController, readout_times, sequences, diffusivity)

Computes the timepoints at which the simulation will be evaluated.

The following timepoints will always be included
- Any multiple of the TR
- Any RF pulse
- Any change in the gradient strength

Additional timepoints will be added to ensure that:
- at any time the timestep is larger than `gradient_precision` times ``D^{-1/3} (G \\gamma)^{-2/3}``.
- between any 2 RF pulses or readouts the timestep is larger than the time between those events divded by the sample_frequency
"""
function get_times(time_controller::TimeController, readout_times::Vector{Float}, sequences::Vector{<:Sequence}, diffusivity::Float)

    timepoints = sorted(readout_times)
    tmax = maximum(timepoints)
    for sequence in sequences
        index = 1
        while time(sequence, index) < tmax
            push!(timepoints, time(sequence, index))
            index += 1
        end

        index = 1
        while index * sequence.TR < tmax
            push!(timepoints, index * sequence.TR)
            index += 1
        end

        append!(timepoints, sequence.gradient.times)
    end
    push!(timepoints, zero(Float))
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
        index_readout = searchsortedfirst(readout_times, mean_time)
        max_timestep_displacement = (
            readout_times[index_readout] -
            readout_times == 1 ? 0 : readout_times[index_readout - 1]
        ) / time_controller.sample_displacement
        for sequence in sequences
            control_points = sequence.gradient.times
            index_control = searchsortedfirst(control_points, mean_time)
            mt = (
                control_points[index_control] -
                control_points[index_control - 1]
            ) / sample_frequency
            if mt < max_timestep_displacement
                max_timestep_displacement = mt
            end
        end

        # sample off-resonance field between RF pulse and subsequence RF pulse or readout
        max_timestep_offresonance = Inf64
        next_readout = readout_times[searchsortedfirst(readout_times, mean_time)]
        for sequence in sequences
            next_rf_pulse = next_pulse(sequence, mean_time)
            if next_rf_pulse <= 1
                continue
            end

            mt = (
                min(next_readout, time(sequence, next_rf_pulse)) - 
                time(sequence, next_rf_pulse - 1)
            )
            if mt < max_timestep_offresonance
                max_timestep_offresonance = mt
            end
        end

        max_timestep = min(max_timestep_gradient, max_timestep_displacement, max_timestep_offresonance)
        @show (max_timestep_gradient, max_timestep_displacement, max_timestep_offresonance)

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


"""
Computes the previous or next RF pulse/readout
"""
function _get_event(current, sequence, readout_times; previous=false)
    correction = previous ? -1 : 0
    index = next_pulse(sequence, current) + correction
    time_pulse = index == 0 ? zero(Float) : time(sequence, index)
    index_readout = searchsortedfirst(readout_times, current) + correction
    time_readout = index_readout == 0 ? zero(Float) : readout_times[index_readout]
    @show (current, time_pulse, time_readout)
    @show previous ? max(time_pulse, time_readout) : min(time_pulse, time_readout)
    return previous ? max(time_pulse, time_readout) : min(time_pulse, time_readout)
end