# defining the sequence
abstract type SequenceComponent end

struct RFPulse{T<:AbstractFloat} <: SequenceComponent
    time :: T
    flip_angle :: T
    cf :: T
    sf :: T
    phase :: T
    cp :: T
    sp :: T
    RFPulse(time, flip_angle, phase) = begin
        f = deg2rad(flip_angle)
        p = deg2rad(phase)
        new{typeof(f)}(time, f, cos(f), sin(f), p, cos(p), sin(p))
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

struct Readout{T} <: SequenceComponent
    time :: T
end
Readout(;time=0.) = Readout(time)

apply(pulse :: Readout, orient :: SpinOrientation) = orient


struct InstantGradient{T} <: SequenceComponent
    qvec :: PosVector{T}
    q_origin :: T
    time :: T
end

InstantGradient(; qvec::AbstractVector=[0., 0., 0.], q_origin=0., time :: Real=0.) = InstantGradient(SVector{3}(qvec), q_origin, time)

function apply(pulse :: InstantGradient, orient :: SpinOrientation, pos::PosVector)
    adjustment = (pos ⋅ pulse.qvec) + pulse.q_origin
    SpinOrientation(orient.longitudinal, orient.transverse, orient.phase + adjustment)
end

apply(pulse :: SequenceComponent, orient :: SpinOrientation, pos::PosVector) = apply(pulse, orient)
apply(pulse :: SequenceComponent, spin :: Spin) = Spin(spin.position, apply(pulse, spin.orientation, spin.position), spin.rng)
apply(pulse :: Nothing, orient :: SpinOrientation) = orient
apply(pulse :: SequenceComponent, snap :: Snapshot) = Snapshot(apply.(pulse, snap.spins), span.time)

apply(pulses :: SVector{N, Union{Nothing, <:SequenceComponent}}, spin :: MultiSpin{N, T}) where {N, T<:AbstractFloat} = MultiSpin{N, T}(
    spin.position,
    SVector{N, SpinOrientation{T}}(SpinOrientation{T}[apply(p, o, spin.position) for (p, o) in zip(pulses, spin.orientations)]),
    spin.rng
)

apply(pulses :: SVector{N, Union{Nothing, <:SequenceComponent}}, snap :: MultiSnapshot{N, T}) where {N, T} = MultiSnapshot{N, T}(
    MultiSpin{N, T}[apply(pulses, s) for s in snap.spins],
    snap.time
)