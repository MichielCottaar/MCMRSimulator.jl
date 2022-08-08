# Defining the microstructure seen by the spin

"""
Describes the spatial distribution of the R1, R2, diffusivity, and off-resonance fields.

Construct using [`field`](@ref).
"""
abstract type Field end
(f::Field)(ms :: AbstractVector{Movement}) = sum(m -> f(m) * m.timestep, ms) / sum(m -> m.timestep, ms)

struct ZeroField <: Field
end

(f::ZeroField)(position :: PosVector) = zero(Float)
(f::ZeroField)(m :: Movement) = zero(Float)

struct ConstantField <: Field
    value :: Float
end

(f::ConstantField)(position :: PosVector) = f.value
(f::ConstantField)(m :: Movement) = f.value

struct GradientField <: Field
    gradient :: PosVector
    offset :: Float
end

(f::GradientField)(position :: PosVector) = position ⋅ f.gradient + f.offset
(f::GradientField)(m :: Movement) = f((m.origin + m.destination) / 2.)

"""
Construct a [`Field`](@ref) object to describe the spatial distribution of some parameter (e.g., R1, R2, off-resonance, diffusivity)

# Constructors
    field()

Sets the parameter to zero everywhere.

    field(value::real)

Sets the parameter to a constant value everywhere.

    field(gradient::PosVector, value::real)

Parameter is described by a spatial gradient.
"""
field() = ZeroField()
field(value :: Real) = iszero(value) ? field() : ConstantField(Float(value))
field(gradient :: AbstractVector{<:Real}, value :: Real) = begin
    if all(gradient .== 0) 
        return field(value) 
    else
        gradient = PosVector(gradient)
        value = Float(value)
        return GradientField(gradient, value)
    end
end
field(field :: Field) = field
field(t::Tuple) = field(t...)

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
struct Microstructure{F1 <: Field, F2 <: Field, F3 <: Field, F4 <: Field, G <: Obstructions}
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
        geometry = SVector{length(geometry)}(geometry)
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
