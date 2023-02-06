# defining the sequence
abstract type InstantComponent end

"""
    InstantRFPulse(;time=0., flip_angle=0., phase=0.)

Instantaneous radio-frequency pulse that flips the spins by `flip_angle` degrees in a plane perpendicular to an axis in the x-y plane with an angle of `phase` degrees with respect to the x-axis at given `time`.
Angles are in degrees and the `time` is in milliseconds.
Angles can be retrieved using [`flip_angle`](@ref) and [`phase`](@ref).
Time can be retrieved using [`get_time`](@ref).
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

function Base.show(io::IO, pulse::InstantRFPulse)
    print(io, "InstantRFPulse: t=$(start_time(pulse))ms, θ=$(flip_angle(pulse))°, ϕ=$(phase(pulse))°;")
end

function InstantRFPulse(time, flip_angle, phase)
    f = Float(flip_angle)
    p = Float(phase)
    frad = deg2rad(f)
    prad = deg2rad(p)
    InstantRFPulse(Float(time), f, cos(frad), sin(frad), p, cos(prad), sin(prad))
end

InstantRFPulse(; time=0, flip_angle=0, phase=0) = InstantRFPulse(time, flip_angle, phase)

"""
    phase(instant_pulse)

Returns the angle in the x-y plane of the axis around with the [`InstantRFPulse`](@ref) rotates in degrees.
"""
phase(pulse :: InstantRFPulse) = pulse.phase

"""
    flip_angle(instant_pulse)

Returns the flip angle of the [`InstantRFPulse`](@ref) in degrees.
"""
flip_angle(pulse :: InstantRFPulse) = pulse.flip_angle

"""
    apply!(sequence_component, spin_orientation[, position])

Applies given sequence component to the spin orientation.
This updates the existing spin orientation.
Some pulses (e.g., [`InstantGradient`](@ref)) require positional information as well.

    apply!(sequence_components, spin)
    apply!(sequence_components, snapshot)

Apply all sequence components to the spin orientation in the [`Spin`](@ref) or to all the spins in [`Snapshot`](@ref).
Sequence components (see [`Sequence`](@ref)) can be `nothing` if there is no sequence component at this time.
"""
function apply!(pulse :: InstantRFPulse, spin :: SpinOrientation)
    Bx_init = spin.transverse * cosd(spin.phase)
    By_init = spin.transverse * sind(spin.phase)
    Bxy_parallel  = pulse.cp * Bx_init + pulse.sp * By_init
    Bxy_perp_init = pulse.cp * By_init - pulse.sp * Bx_init

    Bxy_perp = Bxy_perp_init * pulse.cf + spin.longitudinal * pulse.sf
    Bx = Bxy_parallel * pulse.cp - Bxy_perp * pulse.sp
    By = Bxy_perp * pulse.cp + Bxy_parallel * pulse.sp
    spin.longitudinal = spin.longitudinal * pulse.cf - Bxy_perp_init * pulse.sf
    spin.transverse = sqrt(Bx * Bx + By * By)
    spin.phase = rad2deg(atan(By, Bx))
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


"""
    InstantGradient(; qvec=[0, 0, 0], q_origin=0, time=0.)

Infinitely short gradient pulse that encodes phase information given by `qvec` (units: number of rotations/um) and `q_origin` (units: number of rotations).

The number of time a spins at given `position` is rotated is given by `qvec ⋅ position + q_origin`.

The pulse is applied at given `time` (in milliseconds). Retrieve this time using [`get_time`](@ref).
"""
struct InstantGradient <: InstantComponent
    qvec :: PosVector
    q_origin :: Float
    time :: Float
    InstantGradient(qvec, q_origin, time) = new(PosVector(qvec), Float(q_origin), Float(time))
end

function Base.show(io::IO, pulse::InstantGradient)
    print(io, "InstantGradient: t=$(start_time(pulse))ms, q=$(qvec(pulse))/um;")
end

InstantGradient(; qvec::AbstractVector=[0., 0., 0.], q_origin=0., time :: Real=0.) = InstantGradient(SVector{3}(qvec), q_origin, time)
qvec(pulse::InstantGradient) = pulse.qvec
q_origin(pulse::InstantGradient) = pulse.q_origin
qval(pulse::InstantGradient) = norm(qvec(pulse))

function apply!(pulse :: InstantGradient, orient :: SpinOrientation, pos::PosVector)
    adjustment = (pos ⋅ pulse.qvec) + pulse.q_origin
    orient.phase += 360 * adjustment
end

apply!(pulse :: InstantComponent, orient :: SpinOrientation, pos::PosVector) = apply!(pulse, orient)
apply!(pulse :: InstantComponent, spin :: Spin{1}) = apply!(pulse, spin.orientations[1], spin.position)
apply!(pulse :: Nothing, orient :: SpinOrientation) = nothing
apply!(pulse :: Nothing, orient :: SpinOrientation, pos :: PosVector) = nothing

function apply!(pulses :: SVector{N, <:Any}, spin :: Spin{N}) where {N} 
    for index in 1:N
        apply!(
            pulses[index],
            spin.orientations[index],
            spin.position
        )
    end
end

function apply!(pulses :: SVector{N, <:Any}, spins :: AbstractVector{Spin{N}}) where {N}
    for spin in spins
        apply!(pulses, spin)
    end
end

for cls in (:InstantRFPulse, :Readout, :InstantGradient)
    @eval get_time(p::$cls) = p.time
end