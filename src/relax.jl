"""
    relax(spin_orientation, environment, timestep::Real, gradient_off_resonance=0., B0=3.)

Relaxes the [`SpinOrientation`](@ref) within the R1, R2, and off-resonance given by the [`LocalEnvironment`](@ref) over given `timestep`
"""
function relax(orientation :: SpinOrientation, environment :: LocalEnvironment, timestep :: Real, gradient_off_resonance=0., B0=3.)
    @assert timestep >= 0
    if iszero(timestep)
        return orientation
    end
    timestep = Float(timestep)
    B0 = Float(B0)
    SpinOrientation(
        (1 - (1 - orientation.longitudinal) * exp(-environment.R1 * timestep)),
        orientation.transverse * exp(-environment.R2 * timestep),
        (environment.off_resonance * B0 + gradient_off_resonance) * timestep * gyromagnetic_ratio + orientation.phase
    )
end


"""
    transfer(spin1, spin2, fraction)
    transfer(orientation, obstruction)

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


function transfer(orientation :: SpinOrientation, fraction::Float) where {N}
    inv_fraction = 1 - fraction
    SpinOrientation((1 - (1 - orientation.longitudinal) * inv_fraction), orientation.transverse * inv_fraction, orientation.phase)
end

function transfer(orientations :: SVector{N, SpinOrientation}, obstruction :: ObstructionProperties) where {N}
    fraction = obstruction.MT_fraction
    map(o->transfer(o, fraction), orientations)
end