module MRIBuilderToMCMRSimulator

import MCMRSimulator.SequenceParts: gradient_waveform, get_pulses, get_instants, readout_times, ConstantPulse, SequenceEvent, GradientEvent, PulseEvent

import MRIBuilder as mb


function gradient_waveform(seq::mb.Sequence, dim::Int)
    times = Float64[]
    amplitudes = Float64[]
    for (tstart, block) in mb.iter_blocks(seq)
        for (time, full_ampl) in mb.waveform(block)
            push!(times, time + tstart)
            push!(amplitudes, full_ampl[dim])
        end
    end
    return times, amplitudes
end

function get_pulses(seq::mb.Sequence)
    pulses = Tuple{Float64, Float64, Vector{ConstantPulse}}[]
    for (tstart, block) in mb.iter_blocks(seq)
        for key in keys(block)
            component = block[key]
            if component isa Tuple{<:Number, <:mb.Components.RFPulseComponent}
                delay, event = component
                duration = mb.variables.duration(event)
                if duration > 0
                    gen_event = mb.make_generic(event)
                    times = gen_event.time
                    tstep = (times[end] - times[1]) / (length(times) - 1)
                    @assert all(diff(times) .≈ tstep) "Constant step size required for MCMRSimulator"
                    frequency = (gen_event.phase[2:end] - gen_event.phase[1:end-1]) / (tstep * 360)
                    ampl = (gen_event.amplitude[1:end-1] + gen_event.amplitude[2:end]) / 2
                    push!(pulses, (
                        tstart + delay,
                        tstart + delay + duration,
                        [
                            ConstantPulse(a, p, f) for (a, p, f) in 
                            zip(ampl, gen_event.phase[1:end-1], frequency)
                        ]
                    ))
                end
            end
        end
    end
    return pulses
end


function get_instants(seq::mb.Sequence)
    instants = Tuple{Float64, SequenceEvent}[]
    for (tstart, block) in mb.iter_blocks(seq)
        for key in keys(block)
            component = block[key]
            if component isa Tuple{<:Number, <:mb.Components.RFPulseComponent}
                delay, event = component
                if event isa mb.InstantGradient
                    push!(instants, (tstart + delay, GradientEvent(event.qvec)))
                elseif event isa mb.InstantPulse
                    push!(instants, (tstart + delay, PulseEvent(event.flip_angle, event.phase)))
                end
            end
        end
    end
    return instants
end

readout_times(seq::mb.Sequence) = mb.variables.readout_times(seq)

end