"""
    relax!(spin_orientation, timestep, R1, R2, off_resonance)

Updates [`SpinOrientation`] after evolving for `timestep` with given `R1` (1/ms), `R2` (1/ms), and off_resonance (kHz).
"""
function relax!(orientation :: SpinOrientation, timestep :: Real, R1, R2, off_resonance)
    @assert timestep >= 0
    if iszero(timestep)
        return
    end
    timestep = Float(timestep)
    orientation.longitudinal = (1 - (1 - orientation.longitudinal) * exp(-R1 * timestep))
    orientation.transverse *= exp(-R2 * timestep)
    orientation.phase += 360. * off_resonance * timestep
end


"""
    transfer!(spin1, spin2, fraction)
    transfer!(orientation, obstruction)

Transfers `fraction` of the total spin orientation between `spin1` and `spin2`.
Alternatively transfers spin between spin `orientation` and `obstruction`.
"""
function transfer(spin1 :: SpinOrientation, spin2 :: SpinOrientation, fraction :: Real)
    if iszero(fraction)
        return (spin1, spin2)
    end
    inv_fraction = 1 - fraction
    return (
        SpinOrientation(
            spin1.longitudinal * fraction + spin2.longitudinal * inv_fraction,
            spin1.transverse * fraction + spin2.transverse * inv_fraction,
            spin1.phase
        ),
        SpinOrientation(
            spin2.longitudinal * fraction + spin1.longitudinal * inv_fraction,
            spin2.transverse * fraction + spin1.transverse * inv_fraction,
            spin2.phase
        )
    )
end

function transfer(spin1 :: SVector{N, SpinOrientation}, spin2 :: SVector{N, SpinOrientation}, fraction :: Real) where {N}
    if iszero(fraction)
        return (spin1, spin2)
    end
    map(transfer, spin1, spin2)
end

function transfer(spin1 :: Spin{N}, spin2 :: Spin{N}, fraction :: Real) where {N}
    (orient1, orient2) = transfer(spin1.orientations, spin2.orientations, fraction)
    return (
        Spin(spin1.position, orient1, spin1.rng),
        Spin(spin2.position, orient2, spin2.rng),
    )
end


function transfer!(orientation :: SpinOrientation, fraction::Float)
    inv_fraction = 1 - fraction
    orientation.longitudinal = 1 - (1 - orientation.longitudinal) * inv_fraction
    orientation.transverse *= inv_fraction
end

function transfer!(orientations :: SVector{N, SpinOrientation}, obstruction :: ObstructionProperties, timestep :: Float) where {N}
    fraction = 1 - (1 - obstruction.MT_fraction) ^ sqrt(timestep)
    for o in orientations
        transfer!(o, fraction)
    end
end