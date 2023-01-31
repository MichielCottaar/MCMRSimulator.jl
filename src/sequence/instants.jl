# defining the sequence
abstract type InstantComponent end

"""
    InstantRFPulse(;time=0., flip_angle=0., phase=0.)

Instantaneous radio-frequency pulse that flips the spins by `flip_angle` degrees in a plane perpendicular to an axis in the x-y plane with an angle of `phase` degrees with respect to the x-axis at given `time`.
Angles are in degrees (stored internally as radians) and the `time` is in milliseconds.
"""
struct InstantRFPulse <: InstantComponent
    time :: Float
    flip_angle :: Float
    cf :: Float
    sf :: Float
    phase :: Float
    cp :: Float
    sp :: Float
end

function InstantRFPulse(time, flip_angle, phase)
    f = Float(deg2rad(flip_angle))
    p = Float(deg2rad(phase))
    InstantRFPulse(Float(time), f, cos(f), sin(f), p, cos(p), sin(p))
end

InstantRFPulse(; time=0, flip_angle=0, phase=0) = InstantRFPulse(time, flip_angle, phase)

"""
    phase(pulse)

Returns the angle in the x-y plane of the axis around with the RF pulse rotates in degrees
"""
phase(pulse :: InstantRFPulse) = rad2deg(pulse.phase)
"""
    flip_angle(pulse)

Returns the flip angle of the RF pulse in degrees
"""
flip_angle(pulse :: InstantRFPulse) = rad2deg(pulse.flip_angle)

"""
    apply(sequence_component, spin_orientation[, position])

Applies given sequence component to the spin orientation.
Returns a new [`SpinOrientation`](@ref).
Some pulses (e.g., [`InstantGradient`](@ref)) require positional information as well.

    apply(sequence_components, spin)
    apply(sequence_components, snapshot)

Apply all sequence components to the spin orientation in the [`Spin`](@ref) or to all the spins in [`Snapshot`](@ref) (the return type matches this input type).
Sequence components (see [`Sequence`](@ref)) can be `nothing` if there is no sequence component at this time.
"""
function apply(pulse :: InstantRFPulse, spin :: SpinOrientation)
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

"""
    Readout(;time=0)

Readout the spins at the given `time` (in milliseconds) each TR
"""
struct Readout
    time :: Float
    Readout(time) = new(Float(time))
end
Readout(;time=0.) = Readout(time)

apply(pulse :: Readout, orient :: SpinOrientation) = orient


"""
    InstantGradient(; qvec=[0, 0, 0], q_origin=0, time=0.)

Infinitely short gradient pulse that encodes phase information given by `qvec` and `q_origin`.

The phase offset at every `position` is given by `qvec ⋅ position + q_origin`.

The pulse is applied at given `time` (in milliseconds). every TR.
"""
struct InstantGradient <: InstantComponent
    qvec :: PosVector
    q_origin :: Float
    time :: Float
    InstantGradient(qvec, q_origin, time) = new(PosVector(qvec), Float(q_origin), Float(time))
end

InstantGradient(; qvec::AbstractVector=[0., 0., 0.], q_origin=0., time :: Real=0.) = InstantGradient(SVector{3}(qvec), q_origin, time)
qval(pulse::InstantGradient) = norm(pulse.qvec)

function apply(pulse :: InstantGradient, orient :: SpinOrientation, pos::PosVector)
    adjustment = (pos ⋅ pulse.qvec) + pulse.q_origin
    SpinOrientation(orient.longitudinal, orient.transverse, orient.phase + adjustment)
end

apply(pulse :: InstantComponent, orient :: SpinOrientation, pos::PosVector) = apply(pulse, orient)
apply(pulse :: InstantComponent, spin :: Spin{1}) = Spin(spin.position, SVector{1}([apply(pulse, spin.orientations[1], spin.position)]), spin.rng)
apply(pulse :: Nothing, orient :: SpinOrientation) = orient
apply(pulse :: Nothing, orient :: SpinOrientation, pos :: PosVector) = orient
apply(pulse :: InstantComponent, snap :: Snapshot{1}) = Snapshot(apply.(pulse, snap.spins), span.time)

apply(pulses :: SVector{N, <:Any}, spin :: Spin{N}) where {N} = Spin{N}(
    spin.position,
    SVector{N, SpinOrientation}(SpinOrientation[apply(p, o, spin.position) for (p, o) in zip(pulses, spin.orientations)]),
    spin.rng
)

function apply(pulses :: SVector{N, <:Any}, spins :: AbstractVector{Spin{N}}) where {N}
    map(s->apply(pulses, s), spins)
end

apply(pulses :: SVector{N, <:Any}, snap :: Snapshot{N}) where {N} = Snapshot{N}(
    apply(pulses, snap.spins),
    snap.time
)

for cls in (:InstantRFPulse, :Readout, :InstantGradient)
    @eval get_time(p::$cls) = p.time
end