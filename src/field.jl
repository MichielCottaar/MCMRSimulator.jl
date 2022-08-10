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