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

apply(pulse :: RFPulse, spin :: Spin) = Spin(spin.position, apply_pulse(pulse, spin.orientation))

struct Sequence
    pulses :: Vector{SequenceComponent}
    TR :: Real
    function Sequence(pulses::Vector{T}, TR :: Real) where T <: SequenceComponent
        result = new(sort(pulses, by=x->x.time), TR)
        if length(result.pulses) > 0
            @assert result.pulses[end].time <= TR
        end
        result
    end
end

Sequence(TR :: Real) = Sequence(SequenceComponent[], TR)
Base.getindex(s :: Sequence, index :: Integer) = s.pulses[((index - 1) % length(s.pulses)) + 1]

function time(sequence :: Sequence, index :: Integer)
    if length(sequence.pulses) == 0
        return Inf
    end
    nTR = (index - 1) ÷ length(sequence.pulses)
    return nTR * sequence.TR + sequence[index].time
end
