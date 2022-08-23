"The off-resonance, R2, and R1 values at a single point in space"
struct LocalEnvironment
    off_resonance :: Float
    R2 :: Float
    R1 :: Float
end

"""
    Microstructure(; off_resonance=0., R2=0., R1=0., diffusivity=0., geometry=Obstruction[])

Describes the microstructure of the tissue.

This describes both the spatial distribution of the `R1`, `R2`, `off-resonance`, and `diffusivity` parameters (as [`Field`](@ref) objects)
as well as the spatial `geometry` constraining the diffusion (as a sequence of [`Obstruction`](@ref) objects).

The `R1`, `R2`, `off-resonance`, and `diffusivity` parameters can be defined as one of:
- `value::Number`: constant value across the microstructure.
- `(gradient::PosVector, value::Number)`: gradient across the microstructure.
- `field::Field`: as generated using [`field`](@ref).
"""
struct Microstructure{F1 <: Field, F2 <: Field, F3 <: Field, F4 <: Field, G <: Tuple}
    off_resonance :: F1  # in ppm
    R2 :: F2
    R1 :: F3
    diffusivity :: F4
    geometry :: G
    function Microstructure(;off_resonance=0., R2=0., R1=0., diffusivity=0., geometry=Obstruction[]) 
        if isa(geometry, Obstruction)
            geometry = [geometry]
        end
        R1 = field(R1)
        R2 = field(R2)
        diffusivity = field(diffusivity)
        off_resonance = field(off_resonance)
        if isa(geometry, Obstruction)
            geometry = tuple(geometry)
        else
            geometry = tuple(geometry...)
        end
        if isa(diffusivity, ZeroField) && length(geometry) > 0
            @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
        end
        new{typeof(off_resonance), typeof(R2), typeof(R1), typeof(diffusivity), typeof(geometry)}(
            off_resonance, R2, R1, diffusivity, geometry)
    end
end

(micro::Microstructure)(position) = LocalEnvironment(
    micro.off_resonance(position),
    micro.R2(position),
    micro.R1(position),
)

"""
    relax(orientation::SpinOrientation, environment::LocalEnvironment, timestep::Real, B0=3.)

Relaxes the spin `orientation` within the R1, R2, and off-resonance given by the local `environment` over time `timestep`
"""
function relax(orientation :: SpinOrientation, environment :: LocalEnvironment, timestep :: Real, B0=3.)
    @assert timestep >= 0
    if iszero(timestep)
        return orientation
    end
    timestep = Float(timestep)
    B0 = Float(B0)
    SpinOrientation(
        (1 - (1 - orientation.longitudinal) * exp(-environment.R1 * timestep)),
        orientation.transverse * exp(-environment.R2 * timestep),
        environment.off_resonance * timestep * gyromagnetic_ratio * B0 + orientation.phase
    )
end
