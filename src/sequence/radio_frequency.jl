"""
    RFPulse(times, amplitudes[, phases]; phase0=0., off_resonance=0)

Creates a radio-frequency block between times `t0` and `t1`.
At time `t0` the phase of the transverse B1 magnetic field is equal to `phase0` (given in degrees).

The amplitude and frequency modulation profile of the RF pulse are provided as functions:
- The amplitude represents the strength of the magnetic field in rad/ms the transverse plane as a function of time (between t0 and t1). This should always be positive.
- The frequency represents the offset from the Larmor frequency in rad/ms as a function of time (between t0 and t1). This can be negative or posituT as a function of time (between t0 and t1). This can be negative or posituT as a function of time (between t0 and t1). This can be negative or posituT as a function of time (between t0 and t1). This can be negative or posituT as a function of time (between t0 and t1). This can be negative or positive.
"""
struct RFPulse
    amplitude :: ConcreteShape
    phase :: ConcreteShape
    max_amplitude :: Float
end

function RFPulse(times::AbstractVector, amplitudes::AbstractVector, phases=nothing; phase0=0., off_resonance=0.) 
    if isnothing(phases)
        phases = zeros(Float, length(times))
    end
    phases .+= phase0
    phases = deg2rad.(phases)
    phases .+= times .* off_resonance
    off_resonance_step = @. abs((phases[1:end-1] - phases[2:end]) / (times[1:end-1] - times[2:end]))
    off_resonance_step[isnan.(off_resonance_step)] .= 0
    off_resonance_edge = max.([0, off_resonance_step...], [off_resonance_step..., 0])
    return RFPulse(
        ConcreteShape(times, amplitudes),
        ConcreteShape(times, phases),
        sqrt(maximum(@. amplitudes^2 + off_resonance_edge^2))
    )
end

start_time(pulse::RFPulse) = min(start_time(pulse.amplitude), start_time(pulse.phase))
end_time(pulse::RFPulse) = min(end_time(pulse.amplitude), end_time(pulse.phase))

amplitude(pulse::RFPulse, time::Number) = amplitude(pulse.amplitude, time)
phase(pulse::RFPulse, time::Number) = amplitude(pulse.phase, time)
frequency(pulse::RFPulse, time::Number) = amplitude_derivative(pulse.phase, time)

amplitude(pulse::RFPulse, t1::Number, t2::Number) = amplitude_integral(pulse.amplitude, t1, t2)
phase(pulse::RFPulse, t1::Number, t2::Number) = amplitude_integral(pulse.phase, t1, t2) / (t2 - t1)

function add_TR(pulse::RFPulse, delta_time::Number) 
    RFPulse(
        add_TR(pulse.amplitude, delta_time),
        add_TR(pulse.phase, delta_time),
        pulse.max_amplitude,
    )
end


"""
    effective_pulse(RFPulse, t0, t1)
    effective_pulse(sequence, t0, t1)

Represents the effect of the RF pulse between times `t0` and `t1` as an [`InstantRFPulse`](@ref).
"""
effective_pulse(pulse::RFPulse, t0::Number, t1::Number) = InstantRFPulse(flip_angle=rad2deg(amplitude(pulse, t0, t1)), phase=rad2deg(phase(pulse, t0, t1)))

"""
    constant_pulse(t0, t1, flip_angle; phase0=0., off_resonance=0.)

Creates an RF pulse with a constant amplitude resulting in given `flip_angle` (for spins at `off_resonance`).
"""
function constant_pulse(t0::Number, t1::Number, flip_angle::Number; kwargs...)
    amplitude = Float(deg2rad(flip_angle) / (t1 - t0))
    RFPulse([t0, t0, t1, t1], [0, amplitude, amplitude, 0]; kwargs...)
end