# defining the sequence
abstract type SequenceComponent end

struct RFPulse <: SequenceComponent
    time :: Real
    flip_angle :: Real
    cf :: Real
    sf :: Real
    phase :: Real
    cp :: Real
    sp :: Real
    RFPulse(time, flip_angle, phase) = begin
        f = deg2rad(flip_angle)
        p = deg2rad(phase)
        new(time, f, cos(f), sin(f), p, cos(p), sin(p))
    end
end

RFPulse(; time=0, flip_angle=0, phase=0) = RFPulse(time, flip_angle, phase)

phase(pulse :: RFPulse) = rad2deg(pulse.phase)
flip_angle(pulse :: RFPulse) = rad2deg(pulse.flip_angle)
time(pulse :: RFPulse) = pulse.time

function apply(pulse :: RFPulse, spin :: SpinOrientation)
    Bx_init = spin.transverse * cos(spin.phase)
    By_init = spin.transverse * sin(spin.phase)
    Bxy_parallel  = pulse.cp * Bx_init + pulse.sp * By_init
    Bxy_perp_init = pulse.cp * By_init - pulse.sp * Bx_init

    Bxy_perp = Bxy_perp_init * pulse.cf + spin.longitudinal * pulse.sf
    Bx = Bxy_parallel * pulse.cp - Bxy_perp * pulse.sp
    By = Bxy_perp * pulse.cp + Bxy_parallel * pulse.sp
    SpinOrientation(
        spin.longitudinal * pulse.cf - Bxy_perp_init * pulse.sf,
        sqrt(Bx * Bx + By * By),
        atan(By, Bx)
    )
end

struct Readout <: SequenceComponent
    time :: Real
end
Readout(;time=0.) = Readout(time)

apply(pulse :: Readout, orient :: SpinOrientation) = orient


struct InstantGradient <: SequenceComponent
    qvec :: PosVector
    q_origin :: Real
    time :: Real
end

InstantGradient(; qvec::AbstractVector=[0., 0., 0.], q_origin=0., time :: Real=0.) = InstantGradient(SVector{3}(qvec), q_origin, time)

function apply(pulse :: InstantGradient, orient :: SpinOrientation, pos::PosVector)
    adjustment = (pos ⋅ pulse.qvec) + pulse.q_origin
    SpinOrientation(orient.longitudinal, orient.transverse, orient.phase + adjustment)
end

apply(pulse :: SequenceComponent, orient :: SpinOrientation, pos::PosVector) = apply(pulse, orient)
apply(pulse :: SequenceComponent, spin :: Spin) = Spin(spin.position, apply(pulse, spin.orientation, spin.position), spin.rng)
