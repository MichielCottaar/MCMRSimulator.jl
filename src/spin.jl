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
    position :: SVector
    orientation :: SpinOrientation
end
Spin(;position=zero(SVector{3,Float64}), longitudinal=1., transverse=0., phase=0.) = Spin(SVector{3}(position), SpinOrientation(longitudinal, transverse, deg2rad(phase)))
Base.zero(::Type{Spin}) = Spin()

for param in (:longitudinal, :transverse)
    @eval $param(o :: SpinOrientation) = o.$param
end
phase(o :: SpinOrientation) = norm_angle(rad2deg(o.phase))
vector(o :: SpinOrientation) = SA_F64[
    o.transverse * cos(o.phase),
    o.transverse * sin(o.phase),
    o.longitudinal
]
vector2spin(vector :: SVector{3, T}) where T <: Real = SpinOrientation(
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
const Obstructions = Union{Obstruction, AbstractVector{T}} where T <: Obstruction

struct Movement
    origin :: SVector{3, Real}
    destination :: SVector{3, Real}
    timestep :: Real
end