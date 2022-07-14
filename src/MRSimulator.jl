module MRSimulator
using StaticArrays
using LinearAlgebra

mutable struct Spin
    time :: Real
    position :: SVector{3,Real}
    phase :: Real
    Spin(;time=0., position=zero(SVector{3,Real}), phase=0.) = new(time, position, phase)
end

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
end

struct Microstructure
    off_resonance :: Field
    Microstructure(;off_resonance=ZeroField{Float64}()) = new(off_resonance)
end


(micro::Microstructure)(position :: SVector{3,Real}) = LocalEnvironment(
    micro.off_resonance(position)
)

function relax!(spin :: Spin, micro :: Microstructure, timestep :: Real)
    @assert timestep > 0
    environment = micro(spin.position)
    spin.phase += environment.off_resonance * timestep
end

function evolve!(spin :: Spin, micro :: Microstructure, new_time :: Real)
    if spin.time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time > spin.time
        timestep = new_time - spin.time
        relax!(spin, micro, timestep)
        spin.time = new_time
    end
    spin
end

end