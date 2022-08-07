function norm_angle(angle)
    angle = mod(angle, 360)
    if angle > 180
        angle -= 360
    end
    angle
end

"""Immutable version of the Xoshiro random number generator state

Used to store the current state in the Spin object.
To evolve the spin in a predictable manner set the seed using `copy!(spin.rng, Random.TaskLocalRNG)`.
"""
struct FixedXoshiro
    s0::UInt64
    s1::UInt64
    s2::UInt64
    s3::UInt64
    function FixedXoshiro(state::Random.Xoshiro)
        FixedXoshiro(state.s0, state.s1, state.s2, state.s3)
    end
    FixedXoshiro(s0::Integer, s1::Integer, s2::Integer, s3::Integer) = new(s0, s1, s2, s3)
end

FixedXoshiro(seed=nothing) = FixedXoshiro(Random.Xoshiro(seed))
Random.Xoshiro(rng::FixedXoshiro) = Random.Xoshiro(rng.s0, rng.s2, rng.s2, rng.s3)
function Base.copy!(dst::Random.TaskLocalRNG, src::FixedXoshiro)
    copy!(dst, Random.Xoshiro(src))
    dst
end

struct SpinOrientation{T <: AbstractFloat}
    longitudinal :: T
    transverse :: T
    phase :: T
end

struct Spin{T <: AbstractFloat}
    position :: PosVector{T}
    orientation :: SpinOrientation{T}
    rng :: FixedXoshiro
end

struct MultiSpin{N, T <: AbstractFloat}
    position :: PosVector{T}
    orientations :: SVector{N, SpinOrientation{T}}
    rng :: FixedXoshiro
end

Spin(;position=zero(SVector{3,Float64}), longitudinal=1., transverse=0., phase=0., rng=FixedXoshiro()) = Spin(position, SpinOrientation(longitudinal, transverse, deg2rad(phase)), rng)
MultiSpin(spin::Spin, n_sequences) = MultiSpin(spin.position, SVector{n_sequences}(repeat([spin.orientation], n_sequences)), spin.rng)
Base.zero(::Type{Spin}) = Spin()

for param in (:longitudinal, :transverse)
    @eval $param(o :: SpinOrientation) = o.$param
end
for param in (:*, :/)
    @eval Base.$param(o :: SpinOrientation, number :: Real) = SpinOrientation(
        $param(o.longitudinal, number),
        $param(o.transverse, number),
        o.phase
    )
    @eval Base.$param(s :: Spin, number :: Real) = Spin(s.position, $param(s.orientation, number))
end
phase(o :: SpinOrientation) = norm_angle(rad2deg(o.phase))
vector(o :: SpinOrientation) = SA_F64[
    o.transverse * cos(o.phase),
    o.transverse * sin(o.phase),
    o.longitudinal
]
vector2spin(vector :: AbstractVector) = SpinOrientation(
    vector[3],
    sqrt(vector[1] * vector[1] + vector[2] * vector[2]),
    atan(vector[2], vector[1]),
)

for param in (:longitudinal, :transverse, :phase, :vector)
    @eval $param(s :: Spin) = $param(s.orientation)
end

position(s :: Spin) = s.position

struct Snapshot
    spins :: Vector{Spin}
    time :: Real
end
Snapshot(nspins :: Int, time :: Real) = Snapshot(zeros(Spin, nspins), time)
time(s :: Snapshot) = s.time

function vector(s :: Snapshot)
    sum(vector.(s.spins))
end
for param in (:longitudinal, :transverse, :phase)
    @eval $param(s :: Snapshot) = $param(vector2spin(vector(s)))
end
Base.getindex(s::Snapshot, i::Int) = s.spins[i]
Base.length(s::Snapshot) = length(s.spins)
Base.iterate(s::Snapshot) = iterate(s.spins)
Base.iterate(s::Snapshot, state) = iterate(s.spins, state)

abstract type Obstruction end
const Obstructions{N, T} = SVector{N, T} where T <: Obstruction

struct Movement{T<:AbstractFloat}
    origin :: PosVector{T}
    destination :: PosVector{T}
    timestep :: T
end