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


struct LocalEnvironment
    off_resonance :: Real
    R2 :: Real
    R1 :: Real
end

struct Microstructure
    off_resonance :: Field  # in ppm
    R2 :: Field
    R1 :: Field
    Microstructure(;off_resonance=ZeroField{Float64}(), R2=ZeroField{Float64}(), R1=ZeroField{Float64}()) = new(off_resonance, R2, R1)
end

(micro::Microstructure)(position :: SVector{3,Real}) = LocalEnvironment(
    micro.off_resonance(position),
    micro.R2(position),
    micro.R1(position),
)

field() = ZeroField()
field(value :: Real) = real == zero(typeof(value)) ? ZeroField() : ConstantField(value)
field(gradient :: AbstractVector, value :: Real) = begin
    @assert size(gradient) == (3,)
    all(gradient .== 0) ? field(value) : GradientField(gradient, value)
end

function relax(orient :: SpinOrientation, env :: LocalEnvironment, timestep :: Real, B0=3.)
    @assert timestep > 0
    SpinOrientation(
        (1 - (1 - orient.longitudinal) * exp(-env.R1 * timestep)),
        orient.transverse * exp(-env.R2 * timestep),
        env.off_resonance * timestep * gyromagnetic_ratio * B0 + orient.phase
    )

end
