"""
    Microstructure(; off_resonance=0., R2=0., R1=0., diffusivity=0., geometry=Obstruction[])

Describes the microstructure of the tissue.

This describes both the spatial distribution of the `R1` (in 1/ms), `R2` (in 1/ms), and `off-resonance` (in kHz) parameters (as [`Field`](@ref) objects),
the diffusivity (in um^2/ms) as a number,
as well as the spatial `geometry` constraining the diffusion (as a sequence of [`Obstruction`](@ref) objects).

The `R1`, `R2`, and `off-resonance` parameters can be defined as one of:
- `value::Number`: constant value across the microstructure.
- `(gradient::PosVector, value::Number)`: gradient across the microstructure.
- `field::Field`: as generated using [`field`](@ref).
"""
struct Microstructure{F1 <: Field, F2 <: Field, F3 <: Field, G <: Tuple}
    off_resonance :: F1
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