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


"""
    transfer!(orientation, MT_fraction)

Loses `MT_fraction` spin from `orientation`.
"""
function transfer!(orientation :: SpinOrientation, fraction::Float)
    inv_fraction = 1 - fraction
    orientation.longitudinal = 1 - (1 - orientation.longitudinal) * inv_fraction
    orientation.transverse *= inv_fraction
end