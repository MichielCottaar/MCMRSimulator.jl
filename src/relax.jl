module Relax

import StaticArrays: SVector
import Rotations
import MRIBuilder: LinearSequence, SequencePart
import ..Constants: gyromagnetic_ratio
import ..Spins: Spin, SpinOrientation, stuck, R1, R2, off_resonance, stuck_to, orientation
import ..Geometries.Internal: susceptibility_off_resonance, MRIProperties
import ..Methods: B0
import ..Properties: GlobalProperties

"""
    relax!(spin_orientation, timestep, R1, R2)

Updates [`SpinOrientation`] after evolving for `timestep` with given `R1` (1/ms) and `R2` (1/ms).
"""
function relax!(orientation :: SpinOrientation, timestep :: Real, R1, R2)
    @assert timestep >= 0
    if iszero(timestep)
        return
    end
    timestep = Float64(timestep)
    orientation.longitudinal = (1 - (1 - orientation.longitudinal) * exp(-R1 * timestep))
    orientation.transverse *= exp(-R2 * timestep)
end


"""
    relax!(spin, sequence_parts, simulation, t1, t2)

Evolves a `spin` during part of a sequence represented by one or more [`SequencePart`](@ref) (one for each [`SpinOrientation`](@ref) in the [`Spin`](@ref)).
The spin is evolved for the duration from `t1` and `t2`, where `t=0` corresponds to the start of the sequence part and `t=1` to the end.
The spin is assumed to be stationary during this period.

The spin will precess around an effective RF pulse and relax with a given `R1` and `R2`.
The R1 and R2 are determined based on the [`Simulation`](@ref) parameters, potentially influenced by the geometry.
The effective RF pulse is determined by the sequence and any changes in the off-resonance field due to the geometry or set in the [`Simulation`](@ref).
"""
function relax!(spin::Spin{N}, parts::SVector{N, SequencePart}, simulation, t1::Number, t2::Number) where {N}
    props = MRIProperties(simulation.geometry, simulation.inside_geometry, simulation.properties, spin.position, spin.reflection)

    off_resonance_unscaled = susceptibility_off_resonance(simulation, spin) * 1e-6 * gyromagnetic_ratio

    relax_time = 1 / max(props.R1, props.R2)
    for (orient, part) in zip(spin.orientations, parts)
        off_resonance_kHz = off_resonance_unscaled * B0(part) + props.off_resonance
        if iszero(part.rf_amplitude)
            relax!(orient, step_size, props.R1, props.R2)
            grad = (spin.position - part.gradient_origin) ⋅ part.gradient(t_start + step_ratio/2)
            orient.phase = orient.phase + (grad + off_resonance_kHz) * 360 * (t2 - t1) * part.duration
        else
            nsteps = isinf(relax_time) ? 1 : Int(div(part.duration * (t2 - t1), relax_time / 10, RoundUp))
            step_ratio = (t2 - t1) / nsteps
            step_size = part.duration * step_ratio
            relax!(orient, step_size/2, props.R1, props.R2)
            t_end = t1
            for index in 1:nsteps
                t_start = t_end
                t_end += step_ratio
                grad = (spin.position - part.gradient_origin) ⋅ part.gradient(t_start + step_ratio/2)

                # rotation step
                rotation = RF_pulse_rotation(part, t_start, t_end, grad + off_resonance_kHz)
                new_orient = SpinOrientation(rotation * orientation(orient))
                orient.longitudinal = new_orient.longitudinal
                orient.transverse = new_orient.transverse
                orient.phase = new_orient.phase

                relax!(orient, index == nsteps ? step_size/2 : step_size, props.R1, props.R2)
            end
        end
    end
end


"""
    transfer!(orientation, MT_fraction)

Loses `MT_fraction` spin from `orientation`.
"""
function transfer!(orientation :: SpinOrientation, fraction::Float64)
    inv_fraction = 1 - fraction
    orientation.longitudinal = 1 - (1 - orientation.longitudinal) * inv_fraction
    orientation.transverse *= inv_fraction
end


function RF_pulse_rotation(part::SequencePart, t1, t2, off_resonance)
    frequency = (part.rf_phase.final - part.rf_phase.start) / (360 * part.duration)  # in kHz
    flip_angle = (part.rf_amplitude(t1) + part.rf_amplitude(t2)) * π * part.duration * (t2 - t1)
    phase = part.rf_phase(t1)
    Rotations.RotationVec(
        flip_angle * cosd(phase),
        flip_angle * sind(phase),
        off_resonance - frequency
    )
end

end