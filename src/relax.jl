"""
    relax!(spin_orientation, timestep, mri, additional_off_resonance=0)

Updates [`SpinOrientation`] after evolving for `timestep` with given `R1(mri)` (1/ms), `R2(mri)` (1/ms), and `off_resonance(mri) + additional_off_resonance` (kHz).
"""
function relax!(orientation :: SpinOrientation, timestep :: Real, mri::MRIProperties, additional_off_resonance::Float=zero(Float))
    @assert timestep >= 0
    if iszero(timestep)
        return
    end
    timestep = Float(timestep)
    orientation.longitudinal = (1 - (1 - orientation.longitudinal) * exp(-R1(mri) * timestep))
    orientation.transverse *= exp(-R2(mri) * timestep)
    orientation.phase += 360. * (off_resonance(mri) + additional_off_resonance) * timestep
end


function relax!(spin::Spin{N}, parts::SVector{N, SequencePart}, geometry::Geometry, t1::Number, t2::Number, props::MRIProperties) where {N}
    if stuck(spin)
        use_mri = surface_MRI_properties(Collision(spin), props)
    else
        use_mri = inside_MRI_properties(geometry, spin.position, props)
    end
    off_resonance_unscaled = off_resonance(geometry, spin.position) * 1e-6 * gyromagnetic_ratio
    tmean = (t1 + t2) / 2
    for (orient, part) in zip(spin.orientations, parts)
        pulse = effective_pulse(part, t1, t2)
        off_resonance_kHz = off_resonance_unscaled * B0(part)

        grad = gradient(spin.position, part, t1, tmean)
        relax!(orient, (t2 - t1) * part.total_time/2, use_mri, grad + off_resonance_kHz)
        apply!(pulse, orient)
        grad = gradient(spin.position, part, tmean, t2)
        relax!(orient, (t2 - t1) * part.total_time/2, use_mri, grad + off_resonance_kHz)
    end
end


"""
    transfer!(orientation, MT_fraction)

Loses `MT_fraction` spin from `orientation`.
"""
function transfer!(orientation :: SpinOrientation, fraction::Float)
    inv_fraction = 1 - fraction
    orientation.longitudinal = 1 - (1 - orientation.longitudinal) * inv_fraction
    orientation.transverse *= inv_fraction
end