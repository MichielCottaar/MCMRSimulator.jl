# Defining the microstructure seen by the spin

abstract type Field{T} end
(f::Field)(ms :: AbstractVector{Movement}) = sum(m -> f(m) * m.timestep, ms) / sum(m -> m.timestep, ms)

struct ZeroField{T} <: Field{T}
end

(f::ZeroField{T})(position :: PosVector) where {T} = zero(T)
(f::ZeroField{T})(m :: Movement) where {T} = zero(T)

struct ConstantField{T} <: Field{T}
    value :: T
end

(f::ConstantField)(position :: PosVector) = f.value
(f::ConstantField)(m :: Movement) = f.value

struct GradientField{T} <: Field{T}
    gradient :: SVector{3,T}
    offset :: T
end

(f::GradientField)(position :: PosVector) = position ⋅ f.gradient + f.offset
(f::GradientField)(m :: Movement) = f((m.origin + m.destination) / 2.)

field() = field(Float64)
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

struct LocalEnvironment{T <: Real}
    off_resonance :: T
    R2 :: T
    R1 :: T
end

struct Microstructure{F1 <: Field, F2 <: Field, F3 <: Field, F4 <: Field, G <: Obstructions}
    off_resonance :: F1  # in ppm
    R2 :: F2
    R1 :: F3
    diffusivity :: F4
    geometry :: G
    function Microstructure(;off_resonance=field(), R2=field(), R1=field(), diffusivity=field(), geometry=Obstruction[]) 
        if isa(geometry, Obstruction)
            geometry = [geometry]
        end
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
