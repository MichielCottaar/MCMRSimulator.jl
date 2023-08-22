module Relax

import StaticArrays: SVector
import ..Constants: gyromagnetic_ratio
import ..Spins: Spin, SpinOrientation, stuck, R1, R2, off_resonance, stuck_to
import ..Sequences: SequencePart, effective_pulse, gradient, apply!
import ..Geometries.Internal: susceptibility_off_resonance, MRIProperties
import ..Methods: B0
import ..Properties: GlobalProperties

"""
    relax!(spin_orientation, timestep, mri, additional_off_resonance=0)

Updates [`SpinOrientation`] after evolving for `timestep` with given `R1` (1/ms), `R2` (1/ms), and `off_resonance` (kHz).
"""
function relax!(orientation :: SpinOrientation, timestep :: Real, R1, R2, off_resonance)
    @assert timestep >= 0
    if iszero(timestep)
        return
    end
    timestep = Float64(timestep)
    orientation.longitudinal = (1 - (1 - orientation.longitudinal) * exp(-R1 * timestep))
    orientation.transverse *= exp(-R2 * timestep)
    orientation.phase += 360. * off_resonance * timestep
end


function relax!(spin::Spin{N}, parts::SVector{N, SequencePart}, simulation, t1::Number, t2::Number) where {N}
    props = MRIProperties(simulation.geometry, simulation.inside_geometry, simulation.properties, spin.position, spin.reflection)

    if stuck(spin)
        off_resonance_unscaled = susceptibility_off_resonance(simulation.susceptibility, spin.position, spin.reflection.inside) * 1e-6 * gyromagnetic_ratio
    else
        off_resonance_unscaled = susceptibility_off_resonance(simulation.susceptibility, spin.position) * 1e-6 * gyromagnetic_ratio
    end
    tmean = (t1 + t2) / 2
    for (orient, part) in zip(spin.orientations, parts)
        pulse = effective_pulse(part, t1, t2)
        off_resonance_kHz = off_resonance_unscaled * B0(part) + props.off_resonance

        grad = gradient(spin.position, part, t1, tmean)
        relax!(orient, (t2 - t1) * part.total_time/2, props.R1, props.R2, grad + off_resonance_kHz)
        apply!(pulse, orient)
        grad = gradient(spin.position, part, tmean, t2)
        relax!(orient, (t2 - t1) * part.total_time/2, props.R1, props.R2, grad + off_resonance_kHz)
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
end