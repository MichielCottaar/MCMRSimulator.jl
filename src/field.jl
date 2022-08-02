# Defining the microstructure seen by the spin

abstract type Field{T} end
(f::Field)(ms :: AbstractVector{Movement}) = sum(m -> f(m) * m.timestep, ms) / sum(m -> m.timestep, ms)

struct ZeroField{T} <: Field{T}
end

(f::ZeroField{T})(position :: SVector) where {T} = zero(T)
(f::ZeroField{T})(m :: Movement) where {T} = zero(T)

struct ConstantField{T} <: Field{T}
    value :: T
end

(f::ConstantField)(position :: SVector) = f.value
(f::ConstantField)(m :: Movement) = f.value

struct GradientField{T} <: Field{T}
    gradient :: SVector{3,T}
    offset :: T
end

(f::GradientField)(position :: SVector) = position ⋅ f.gradient + f.offset
(f::GradientField)(m :: Movement) = f((m.origin + m.destination) / 2.)

field() = field(typeof(0.))
field(::Type{T}) where T <: Real = ZeroField{T}()
field(value :: Real) = iszero(value) ? field(typeof(value)) : ConstantField(value)
field(gradient :: AbstractVector, value :: Real) = begin
    @assert size(gradient) == (3,)
    if all(gradient .== 0) 
        return field(value) 
    else
        new_type = Base.promote_eltype(gradient, value)
        return GradientField{new_type}(SVector{3,new_type}(gradient), new_type(value))
    end
end

struct LocalEnvironment
    off_resonance :: Real
    R2 :: Real
    R1 :: Real
end

struct Microstructure
    off_resonance :: Field  # in ppm
    R2 :: Field
    R1 :: Field
    diffusivity :: Field
    geometry :: Obstructions
    Microstructure(;off_resonance=field(), R2=field(), R1=field(), diffusivity=field(), geometry=Obstruction[]) = new(off_resonance, R2, R1, diffusivity, geometry)
end

(micro::Microstructure)(position) = LocalEnvironment(
    micro.off_resonance(position),
    micro.R2(position),
    micro.R1(position),
)

function relax(orient :: SpinOrientation, env :: LocalEnvironment, timestep :: Real, B0=3.)
    @assert timestep >= 0
    if timestep == 0
        return orient
    end
    SpinOrientation(
        (1 - (1 - orient.longitudinal) * exp(-env.R1 * timestep)),
        orient.transverse * exp(-env.R2 * timestep),
        env.off_resonance * timestep * gyromagnetic_ratio * B0 + orient.phase
    )
end
