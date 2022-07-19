# Defining the microstructure seen by the spin

abstract type Field{T} end

struct ZeroField{T} <: Field{T}
end

(f::ZeroField{T})(position :: SVector{3,Real}) where {T} = zero(T)


struct ConstantField{T} <: Field{T}
    value :: T
end

(f::ConstantField{T})(position :: SVector{3,Real}) where {T} = f.value

struct GradientField{T} <: Field{T}
    gradient :: SVector{3,T}
    offset :: T
end

(f::GradientField{T})(position :: SVector{3,Real}) where {T} = position ⋅ f.gradient + f.offset

field() = field(typeof(0.))
field(::Type{T}) where T <: Real = ZeroField{T}()
field(value :: Real) = real == zero(typeof(value)) ? field(typeof(T)) : ConstantField(value)
field(gradient :: AbstractVector, value :: Real) = begin
    @assert size(gradient) == (3,)
    all(gradient .== 0) ? field(value) : GradientField(gradient, value)
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
    Microstructure(;off_resonance=field(), R2=field(), R1=field()) = new(off_resonance, R2, R1)
end

(micro::Microstructure)(position :: SVector{3,Real}) = LocalEnvironment(
    micro.off_resonance(position),
    micro.R2(position),
    micro.R1(position),
)

function relax(orient :: SpinOrientation, env :: LocalEnvironment, timestep :: Real, B0=3.)
    @assert timestep > 0
    SpinOrientation(
        (1 - (1 - orient.longitudinal) * exp(-env.R1 * timestep)),
        orient.transverse * exp(-env.R2 * timestep),
        env.off_resonance * timestep * gyromagnetic_ratio * B0 + orient.phase
    )

end
