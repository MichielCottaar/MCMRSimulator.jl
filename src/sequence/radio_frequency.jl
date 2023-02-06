"""
    RFPulse(times, amplitudes[, phases]; phase0=0., off_resonance=0)

Creates a radio-frequency block between times `t0` and `t1`.
At time `t0` the phase of the transverse B1 magnetic field is equal to `phase0` (given in degrees).

The amplitude and frequency modulation profile of the RF pulse are provided as functions:
- The amplitude represents the strength of the magnetic field in kHz the transverse plane as a function of time (between t0 and t1). This should always be positive.
- The phase represents the offset from the Larmor frequency in degrees as a function of time (between t0 and t1). If not provided explicitly, the phase at time `t` is determined by the initial phase (in degrees) and the off-resonance (in kHz).

The functions [`amplitude`](@ref), [`phase`](@ref), and [`off_resonance`](@ref) are used to get the value at a specific `t` or averaged between two times.
"""
struct RFPulse
    amplitude :: ConcreteShape
    phase :: ConcreteShape
    max_amplitude :: Float
end

function RFPulse(times::AbstractVector, amplitudes::AbstractVector, phases=nothing; phase0=0., off_resonance=0.) 
    if isnothing(phases)
        phases = (times .- times[1]) .* (360 * off_resonance) .+ phase0
    end
    off_resonance_step = @. abs((phases[1:end-1] - phases[2:end]) / (360 * (times[1:end-1] - times[2:end])))
    off_resonance_step[isnan.(off_resonance_step)] .= 0
    off_resonance_edge = max.([0, off_resonance_step...], [off_resonance_step..., 0])
    return RFPulse(
        ConcreteShape(times, amplitudes),
        ConcreteShape(times, phases),
        sqrt(maximum(@. amplitudes^2 + off_resonance_edge^2))
    )
end

function Base.show(io::IO, pulse::RFPulse)
    print(io, "RFPulse: t=$(start_time(pulse)) to $(end_time(pulse))ms, θ=$(flip_angle(pulse))°, ϕ=$(phase(pulse, start_time(pulse)))°;")
end

flip_angle(pulse::RFPulse) = amplitude_integral(pulse.amplitude, pulse.amplitude.t0, pulse.amplitude.t1) * 360

start_time(pulse::RFPulse) = min(start_time(pulse.amplitude), start_time(pulse.phase))
end_time(pulse::RFPulse) = max(end_time(pulse.amplitude), end_time(pulse.phase))

"""
    amplitude(rf_pulse, t1[, t2])

Computes the amplitude of the [`RFPulse`](@ref) at time `t1` in kHz.
If `t2` is also provided, the average amplitude between times `t1` and `t2` is returned.
"""
amplitude(pulse::RFPulse, time::Number) = amplitude(pulse.amplitude, time)

"""
    phase(rf_pulse, t1[, t2])

Computes the phase of the [`RFPulse`](@ref) at time `t1` in degrees.
If `t2` is also provided, the average phase between times `t1` and `t2` is returned.
"""
phase(pulse::RFPulse, time::Number) = amplitude(pulse.phase, time)

"""
    off_resonance(rf_pulse, t1[, t2])

Computes the off_resonance of the [`RFPulse`](@ref) at time `t1` in kHz.
If `t2` is also provided, the average off_resonance between times `t1` and `t2` is returned.
"""
off_resonance(pulse::RFPulse, time::Number) = amplitude_derivative(pulse.phase, time) / 360

amplitude(pulse::RFPulse, t1::Number, t2::Number) = amplitude_integral(pulse.amplitude, t1, t2) / (t2 - t1)
phase(pulse::RFPulse, t1::Number, t2::Number) = amplitude_integral(pulse.phase, t1, t2) / (t2 - t1)
off_resonance(pulse::RFPulse, t1::Number, t2::Number) = (amplitude(pulse.phase, t2) - amplitude(pulse.phase, t1)) / ((t2 - t1) * 360)

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
effective_pulse(pulse::RFPulse, t0::Number, t1::Number) = InstantRFPulse(flip_angle=amplitude(pulse, t0, t1) * (t1 - t0) * 360., phase=phase(pulse, t0, t1))

"""
    constant_pulse(t0, t1, flip_angle; phase0=0., off_resonance=0.)

Creates an RF pulse with a constant amplitude resulting in given `flip_angle` (for spins at `off_resonance`).
"""
function constant_pulse(t0::Number, t1::Number, flip_angle::Number; kwargs...)
    amplitude = Float(flip_angle / (360. * (t1 - t0)))
    RFPulse([t0, t1], [amplitude, amplitude]; kwargs...)
end