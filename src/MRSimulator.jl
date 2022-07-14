module MRSimulator
using StaticArrays
using LinearAlgebra

const gyromagnetic_ratio = 267.52218744  # (10^6 rad⋅s^−1⋅T^−1)

mutable struct Spin
    time :: Real
    position :: SVector{3,Real}
    longitudinal :: Real
    transverse :: Real
    phase :: Real
    Spin(;time=0., position=zero(SVector{3,Real}), longitudinal=1., transverse=0., phase=0.) = new(time, position, longitudinal, transverse, phase)
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
    off_resonance :: Field  # in ppm
    Microstructure(;off_resonance=ZeroField{Float64}()) = new(off_resonance)
end


(micro::Microstructure)(position :: SVector{3,Real}) = LocalEnvironment(
    micro.off_resonance(position)
)

function relax!(spin :: Spin, micro :: Microstructure, timestep :: Real, B0=3.)
    @assert timestep > 0
    environment = micro(spin.position)
    spin.phase += environment.off_resonance * timestep * gyromagnetic_ratio * B0
end

function evolve!(spin :: Spin, micro :: Microstructure, new_time :: Real, B0=3.)
    if spin.time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time > spin.time
        timestep = new_time - spin.time
        relax!(spin, micro, timestep, B0)
        spin.time = new_time
    end
    spin
end

end