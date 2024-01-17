module Instants
import LinearAlgebra: ⋅, norm
import StaticArrays: SVector
import ...Spins: Spin, SpinOrientation, Snapshot
import ...Methods: get_time, get_rotation, phase
import ..Methods: start_time, end_time

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
    time :: Float64
    flip_angle :: Float64
    cf :: Float64
    sf :: Float64
    phase :: Float64
    cp :: Float64
    sp :: Float64
end

function Base.show(io::IO, pulse::InstantRFPulse)
    print(io, "InstantRFPulse: t=$(start_time(pulse))ms, θ=$(flip_angle(pulse))°, ϕ=$(phase(pulse))°;")
end

function InstantRFPulse(time, flip_angle, phase)
    f = Float64(flip_angle)
    p = Float64(phase)
    frad = deg2rad(f)
    prad = deg2rad(p)
    InstantRFPulse(Float64(time), f, cos(frad), sin(frad), p, cos(prad), sin(prad))
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

    Bxy_perp = Bxy_perp_init * pulse.cf - spin.longitudinal * pulse.sf
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
    time :: Float64
    Readout(time) = new(Float64(time))
end
Readout(;time=0.) = Readout(time)


"""
    InstantGradient(; qvec=[0, 0, 0], q_origin=0, origin=[0, 0, 0], time=0., apply_bvec=false)

Infinitely short gradient pulse that encodes phase information given by `qvec` (units: rad/um) and `q_origin` (units: rad) or `origin` (units: um).

The number of time a spins at given `position` is rotated is given by `qvec ⋅ position + q_origin` or `qvec ⋅ (position - origin)`.

The pulse is applied at given `time` (in milliseconds). Retrieve this time using [`get_time`](@ref).

If `apply_bvec` is set to true, the gradients will be rotated with the user-provided bvecs file (using [`rotate_bvec`](@ref)).
This should be true for diffusion-weighted gradients, but will typically be false for crusher gradients.
"""
struct InstantGradient <: InstantComponent
    qvec :: SVector{3, Float64}
    origin :: SVector{3, Float64}
    time :: Float64
    apply_bvec :: Bool
    InstantGradient(qvec, origin, time, apply_bvec) = new(SVector{3, Float64}(qvec), SVector{3, Float64}(origin), Float64(time), Bool(apply_bvec))
end

function Base.show(io::IO, pulse::InstantGradient)
    print(io, "InstantGradient: t=$(start_time(pulse))ms, q=$(qvec(pulse))rad/um;")
end

function InstantGradient(; qvec::AbstractVector=[0., 0., 0.], q_origin=nothing, origin=nothing, time :: Real=0., apply_bvec=false)
    if isnothing(origin)
        if isnothing(q_origin)
            origin = zero(SVector{3, Float64})
        else
            dist = q_origin ./ norm(qvec)
            origin = dist .* qvec
        end
    elseif ~isnothing(q_origin)
        @assert q_origin ≈ origin ⋅ qvec
    end
    InstantGradient(SVector{3}(qvec), origin, time, apply_bvec)
end
qvec(pulse::InstantGradient) = pulse.qvec
q_origin(pulse::InstantGradient) = pulse.origin ⋅ pulse.qvec
qval(pulse::InstantGradient) = norm(qvec(pulse))

function apply!(pulse :: InstantGradient, orient :: SpinOrientation, pos::SVector{3, Float64})
    adjustment = (pos ⋅ pulse.qvec) + q_origin(pulse)
    orient.phase += 360 * adjustment / 2π
end

apply!(pulse :: InstantComponent, orient :: SpinOrientation, pos::SVector{3, Float64}) = apply!(pulse, orient)
apply!(pulse :: InstantComponent, spin :: Spin{1}) = apply!(pulse, spin.orientations[1], spin.position)
apply!(pulse :: Nothing, orient :: SpinOrientation) = nothing
apply!(pulse :: Nothing, orient :: SpinOrientation, pos :: SVector{3, Float64}) = nothing

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

function rotate_bvec(gradient::InstantGradient, bvec)
    rotation = get_rotation(bvec, 3)
    InstantGradient(
        rotation * gradient.qvec,
        gradient.origin,
        gradient.time,
        gradient.apply_bvec
    )
end

end