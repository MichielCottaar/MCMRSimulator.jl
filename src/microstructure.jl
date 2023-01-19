"The off-resonance, R2, and R1 values at a single point in space"
struct LocalEnvironment
    off_resonance :: Float
    R2 :: Float
    R1 :: Float
end

"""
    Microstructure(; off_resonance=0., R2=0., R1=0., diffusivity=0., geometry=Obstruction[])

Describes the microstructure of the tissue.

This describes both the spatial distribution of the `R1`, `R2`, and `off-resonance` parameters (as [`Field`](@ref) objects),
the diffusivity as a float,
as well as the spatial `geometry` constraining the diffusion (as a sequence of [`Obstruction`](@ref) objects).

The `R1`, `R2`, and `off-resonance` parameters can be defined as one of:
- `value::Number`: constant value across the microstructure.
- `(gradient::PosVector, value::Number)`: gradient across the microstructure.
- `field::Field`: as generated using [`field`](@ref).
"""
struct Microstructure{F1 <: Field, F2 <: Field, F3 <: Field, G <: Tuple}
    off_resonance :: F1  # in ppm
    R2 :: F2
    R1 :: F3
    diffusivity :: Float
    geometry :: G
    function Microstructure(;off_resonance=0., R2=0., R1=0., diffusivity=0., geometry=Obstruction[]) 
        if isa(geometry, Obstruction)
            geometry = [geometry]
        end
        R1 = field(R1)
        R2 = field(R2)
        diffusivity = Float(diffusivity)
        off_resonance = field(off_resonance)
        if isa(geometry, Obstruction)
            geometry = tuple(geometry)
        else
            geometry = tuple(geometry...)
        end
        if iszero(diffusivity) && length(geometry) > 0
            @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
        end
        new{typeof(off_resonance), typeof(R2), typeof(R1), typeof(geometry)}(
            off_resonance, R2, R1, diffusivity, geometry)
    end
end

(micro::Microstructure)(position) = LocalEnvironment(
    micro.off_resonance(position) + off_resonance(micro.geometry, position),
    micro.R2(position),
    micro.R1(position),
)