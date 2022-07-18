# defining the sequence
abstract type SequenceComponent end

struct EndSequence <: SequenceComponent
    time :: Real
end

apply(s :: EndSequence, spin :: Spin) = spin

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

phase(pulse :: RFPulse) = rad2deg(pulse.phase)
flip_angle(pulse :: RFPulse) = rad2deg(pulse.flip_angle)
time(pulse :: RFPulse) = pulse.time

function apply_pulse(pulse :: RFPulse, spin :: SpinOrientation)
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

apply(pulse :: RFPulse, spin :: Spin) = Spin(spin.time, spin.position, apply_pulse(pulse, spin.orientation))

struct Sequence
    pulses :: Vector{SequenceComponent}
    function Sequence(pulses::Vector{SequenceComponent}, end_sequence :: Real)
        new(sort([pulses; EndSequence(end_sequence)], by=x->x.time))
    end
end

Sequence(end_sequence :: Real) = Sequence(SequenceComponent[], end_sequence)
Base.getindex(s :: Sequence, index :: Integer) = s.pulses[index]

function time(sequence :: Sequence, index :: Integer)
    if index > length(sequence.pulses)
        return Inf
    end
    return sequence[index].time
end
