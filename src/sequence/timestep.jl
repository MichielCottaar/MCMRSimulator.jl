"""
    get_times(simulation::Simulation, readout_times)

Computes the timepoints at which the simulation will be evaluated.

The following timepoints will always be included
- Any multiple of the TR
- Any RF pulse
- Any change in the gradient strength

Additional timepoints will be added to ensure that:
- at any time the timestep is larger than `gradient_precision` times ``D^{-1/3} (G \\gamma)^{-2/3}``.
- between any 2 RF pulses or readouts the timestep is larger than the time between those events divded by the sample_frequency
"""
function get_times(sequences::Vector{<:Sequence}, readout_times::Vector{Float}, diffusivity::Float, gradient_precision::Float, sample_frequency::Int)
    timepoints = copy(readout_times)
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
    end
    timepoints = unique(timepoints)
    sort!(timepoints)

    for (t0, t1) in zip(copy(timepoints[1:end-1]), copy(timepoints[2:end]))
        max_timestep_gradient = Inf64
        if !iszero(gradient_precision)
            for sequence in sequences
                grad = norm(get_gradient(sequence, t0, t1)) * gyromagnetic_ratio * Float(1e-3) * sequence.scanner.B0  # frequency gradient in rad / (ms um)
                mt = gradient_precision * (diffusivity * grad * grad) ^ (-1/3)
                if mt < max_timestep_gradient
                    max_timestep_gradient = mt
                end
            end
        end

        max_timestep_sampling = Inf64
        for sequence in sequences
            mt = (
                _get_event((t0 + t1) / 2, sequence, readout_times; previous=false) - 
                _get_event((t0 + t1) / 2, sequence, readout_times; previous=true)
            ) / sample_frequency
            @show mt
            if mt < max_timestep_sampling
                max_timestep_sampling = mt
            end
        end

        max_timestep = min(max_timestep_gradient, max_timestep_sampling)
        @show max_timestep_gradient
        @show max_timestep_sampling

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