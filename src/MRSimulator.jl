module MRSimulator
import StaticArrays: SA_F64, MVector, SVector
using LinearAlgebra
import Base

const gyromagnetic_ratio = 267.52218744  # (10^6 rad⋅s^−1⋅T^−1)

function norm_angle(angle)
    angle = mod(angle, 360)
    if angle > 180
        angle -= 360
    end
    angle
end

struct SpinOrientation
    longitudinal :: Real
    transverse :: Real
    phase :: Real
end

struct Spin
    time :: Real
    position :: SVector{3,Real}
    orientation :: SpinOrientation
end
Spin(;time=0., position=zero(SVector{3,Real}), longitudinal=1., transverse=0., phase=0.) = Spin(time, position, SpinOrientation(longitudinal, transverse, deg2rad(phase)))

for param in (:longitudinal, :transverse)
    @eval $param(o :: SpinOrientation) = o.$param
    @eval $param(s :: Spin) = $param(s.orientation)
end
phase(o :: SpinOrientation) = norm_angle(rad2deg(o.phase))

phase(s :: Spin) = phase(s.orientation)
time(s :: Spin) = s.time
position(s :: Spin) = s.position

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

# defining the sequence
abstract type SequenceComponent end

struct EndSequence <: SequenceComponent
    time :: Real
end

apply(s :: EndSequence, spin :: Spin) = spin

struct RFPulse <: SequenceComponent
    time :: Real
    flip_angle :: Real
    cf :: Real
    sf :: Real
    phase :: Real
    cp :: Real
    sp :: Real
    RFPulse(time, flip_angle, phase) = begin
        f = deg2rad(flip_angle)
        p = deg2rad(phase)
        new(time, f, cos(f), sin(f), p, cos(p), sin(p))
    end
end

phase(pulse :: RFPulse) = rad2deg(pulse.phase)
flip_angle(pulse :: RFPulse) = rad2deg(pulse.flip_angle)
time(pulse :: RFPulse) = pulse.time

function apply_pulse(pulse :: RFPulse, spin :: SpinOrientation)
    Bx_init = spin.transverse * cos(spin.phase)
    By_init = spin.transverse * sin(spin.phase)
    Bxy_parallel  = pulse.cp * Bx_init + pulse.sp * By_init
    Bxy_perp_init = pulse.cp * By_init - pulse.sp * Bx_init

    Bxy_perp = Bxy_perp_init * pulse.cf + spin.longitudinal * pulse.sf
    Bx = Bxy_parallel * pulse.cp - Bxy_perp * pulse.sp
    By = Bxy_perp * pulse.cp + Bxy_parallel * pulse.sp
    SpinOrientation(
        spin.longitudinal * pulse.cf - Bxy_perp_init * pulse.sf,
        sqrt(Bx * Bx + By * By),
        atan(By, Bx)
    )
end

apply(pulse :: RFPulse, spin :: Spin) = Spin(spin.time, spin.position, apply_pulse(pulse, spin.orientation))

struct Sequence
    pulses :: Vector{SequenceComponent}
    function Sequence(pulses::Vector{SequenceComponent}, end_sequence :: Real)
        new(sort([pulses; EndSequence(end_sequence)], by=x->x.time))
    end
end

Sequence(end_sequence :: Real) = Sequence(SequenceComponent[], end_sequence)
Base.getindex(s :: Sequence, index :: Integer) = s.pulses[index]

function time(sequence :: Sequence, index :: Integer)
    if index > length(sequence.pulses)
        return Inf
    end
    return sequence[index].time
end

function relax(orient :: SpinOrientation, env :: LocalEnvironment, timestep :: Real, B0=3.)
    @assert timestep > 0
    SpinOrientation(
        (1 - (1 - orient.longitudinal) * exp(-env.R1 * timestep)),
        orient.transverse * exp(-env.R2 * timestep),
        env.off_resonance * timestep * gyromagnetic_ratio * B0 + orient.phase
    )

end

function evolve_to_time(spin :: Spin, micro :: Microstructure, new_time :: Real, B0=3.)
    if spin.time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time > spin.time
        timestep = new_time - spin.time
        spin = Spin(
            new_time,
            spin.position,
            relax(spin.orientation, micro(spin.position), timestep, B0)
        )
    end
    spin
end

function evolve(spin :: Spin, micro :: Microstructure, sequence :: Sequence; store_every=1., B0=3.)
    sequence_index = 1
    readout_index = 1
    spins = typeof(spin)[]
    times = MVector{2, Float64}([
        time(sequence, sequence_index),
        (readout_index - 1) * store_every
    ])
    while !isinf(times[1])
        next = argmin(times)
        spin = evolve_to_time(spin, micro, times[next], B0)
        if next == 1
            spin = apply(sequence[sequence_index], spin)
            sequence_index += 1
            times[1] = time(sequence, sequence_index)
        elseif next == 2
            push!(spins, spin)
            readout_index += 1
            times[2] = (readout_index - 1) * store_every
        end
    end
    spins
end

function evolve(spins :: Vector{Spin}, micro :: Microstructure, sequence :: Sequence; store_every=1., B0=3.)
    nspins = length(spins)
    result = fill(Spin[], nspins)
    Threads.@threads for i = 1:nspins
        result[i] = evolve(spins[i], micro, sequence, store_every=store_every, B0=B0)
    end
    hcat(result...)
end
end

