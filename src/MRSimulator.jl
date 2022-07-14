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
end

struct Microstructure
    off_resonance :: Field  # in ppm
    Microstructure(;off_resonance=ZeroField{Float64}()) = new(off_resonance)
end

(micro::Microstructure)(position :: SVector{3,Real}) = LocalEnvironment(
    micro.off_resonance(position)
)

# defining the sequence
abstract type SequenceComponent end

struct RFPulse <: SequenceComponent
    time :: Real
    flip_angle :: Real
    cf :: Real
    sf :: Real
    phase :: Real
    cp :: Real
    sp :: Real
    RFPulse(time, flip_angle, phase) = begin
        f = flip_angle / 180. * π
        p = phase / 180. * π
        new(time, f, cos(f), sin(f), p, cos(p), sin(p))
    end
end

function apply!(pulse :: RFPulse, spin :: Spin)
    @assert pulse.time == spin.time
    Bx_init = spin.transverse * cos(spin.phase)
    By_init = spin.transverse * sin(spin.phase)
    Bxy_parallel  = pulse.cp * Bx_init + pulse.sp * By_init
    Bxy_perp_init = pulse.cp * By_init - pulse.sp * Bx_init

    Bxy_perp = Bxy_perp_init * pulse.cf + spin.longitudinal * pulse.sf
    spin.longitudinal = spin.longitudinal * pulse.cf - Bxy_perp_init * pulse.sf
    Bx = Bxy_parallel * pulse.cp - Bxy_perp * pulse.sp
    By = Bxy_perp * pulse.cp + Bxy_parallel * pulse.sp
    spin.transverse = sqrt(Bx * Bx + By * By)
    spin.phase = atan(By, Bx)
end

struct Sequence
    pulses :: Vector{RFPulse}
    function Sequence(pulses::Vector{RFPulse})
        for c in components
            append!(result, (c.time, c))
        end
        new(sort(pulses, by=x->x.time))
    end
end

function next_time(current_time :: Real, sequence :: Sequence)
    index = searchsortedfirst(sequence.components, current_time, by=x->x.time)
    return sequence[index].time
end

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