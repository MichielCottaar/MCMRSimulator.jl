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

"""
    SpinOrientation(longitudinal, transverse, phase)

The spin orientation. Usually created as part of a [`Spin`](@ref) or [`Multispin`](@ref) object.

This information can be extracted using:
- [`longitudinal`](@ref) to get the spin in the z-direction (equilibrium of 1)
- [`transverse`](@ref) to get the spin in the x-y-plane
- [`phase`](@ref) to get the spin angle in x-y plane (in degrees)
- [`vector`](@ref) to get the length-3 vector
"""
struct SpinOrientation{T <: AbstractFloat}
    longitudinal :: T
    transverse :: T
    phase :: T
end

"""
    Spin(;position=[0, 0, 0], longitudinal=1., transverse=0., phase=0.)

Spin particle with a position and spin orientation (stored as [`SpinOrientation`](@ref)).

A random number generator is stored in the `Spin` object as well, which will be used for evolving the spin into the future in a reproducible manner.

Orientational information can be extracted using:
- [`longitudinal`](@ref) to get the spin in the z-direction (equilibrium of 1)
- [`transverse`](@ref) to get the spin in the x-y-plane
- [`phase`](@ref) to get the spin angle in x-y plane (in degrees)
- [`vector`](@ref) to get a length-3 vector with the spin orientation
- [`position`](@ref) to get a length-3 vector with spin location
"""
struct Spin{T <: AbstractFloat}
    position :: PosVector{T}
    orientation :: SpinOrientation{T}
    rng :: FixedXoshiro
end

Spin(;position=zero(SVector{3,Float64}), longitudinal=1., transverse=0., phase=0., rng=FixedXoshiro()) = Spin(SVector{3}(position), SpinOrientation(longitudinal, transverse, deg2rad(phase)), rng)
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Spin}) = Spin(position=SVector{3}(rand(rng, 3)))
Base.zero(::Type{Spin}) = Spin()

"""
    longitudinal(spin)
    longitudinal(snapshot)

Returns the longitudinal spin (i.e., magnitude aligned with the mangetic field) for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available (i.e., [`MultiSpin`](@ref) or [`MultiSnapshot`](@ref)) an array of longitudinal values is returned with a value for each sequence.
"""
function longitudinal end

"""
    transverse(spin)
    transverse(snapshot)

Returns the transverse spin (i.e., magnitude in the plane perpendicular to the magnetic field) for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available (i.e., [`MultiSpin`](@ref) or [`MultiSnapshot`](@ref)) an array of transverse values is returned with a value for each sequence.
"""
function transverse end

"""
    phase(spin)
    phase(snapshot)

Returns the phase in the x-y plane of the spin for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available (i.e., [`MultiSpin`](@ref) or [`MultiSnapshot`](@ref)) an array of phase values is returned with a value for each sequence.
"""
function phase end

"""
    vector(spin)
    vector(snapshot)

Returns the spin orientation as a length-3 vector for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available (i.e., [`MultiSpin`](@ref) or [`MultiSnapshot`](@ref)) an array of vectors is returned with a value for each sequence.
"""
function vector end

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

"""
Spin particle with a position and a spin orientation stored for every sequence.

This is the equivalent of a [`Spin`](@ref) object used when multiple spins are available.

Orientational information can be extracted using:
- [`longitudinal`](@ref) to get the spin sizes in the z-direction (equilibrium of 1)
- [`transverse`](@ref) to get the spin sizes in the x-y-plane
- [`phase`](@ref) to get the spin angles in x-y plane (in degrees)
- [`vector`](@ref) to get length-3 vectors with the spin orientation
- [`position`](@ref) to get a length-3 vector with spin location
"""
struct MultiSpin{N, T <: AbstractFloat}
    position :: PosVector{T}
    orientations :: SVector{N, SpinOrientation{T}}
    rng :: FixedXoshiro
end

MultiSpin(spin::Spin, n_sequences::Integer) = MultiSpin(spin.position, SVector{n_sequences}(repeat([spin.orientation], n_sequences)), spin.rng)
MultiSpin(n_sequences::Integer; kwargs...) = MultiSpin(Spin(;kwargs...), n_sequences)
"""
    get_sequence(spin, sequence_index)
    get_sequence(snapshot, sequence_index)

Extracts the spin orientation corresponding to a specific sequence, where the sequence index uses the order in which the sequences where provided in the [`Simulation`](@ref).
This converts [`MultiSpin`](@ref) objects into [`Spin`](@ref) objects or [`MultiSnapshot`](@ref) objects into [`Snapshot`](@ref) objects.
"""
get_sequence(spin::MultiSpin, index) = Spin(spin.position, spin.orientations[index], spin.rng)

for param in (:longitudinal, :transverse, :phase, :vector)
    @eval $param(s :: MultiSpin) = $param.(s.orientation)
end


"Returns the position of the spin particle as a vector of length 3."
position(s :: Spin) = s.position
position(s :: MultiSpin) = s.position

"""
Represents the positions and orientations of multiple [`Spin`](@ref) objects at a specific `time`.

This object only keeps track of a single orientation per spin. 
When simulating multiple sequences at once, the [`MultiSnapshot`](@ref) object can be used to keep track of the orientation for each sequene.

Note that times are in milliseconds and positions in micrometer. 
The equilibrium longitudinal spin (after T1 relaxation) is always 1.

# Useful constructors
    Snapshot(positions; time=0., longitudinal=1., transverse=0., phase=0.)
    Snapshot(nspins::Integer; time=0., longitudinal=1., transverse=0., phase=0.)

Creates a new Snapshot at the given `time` with perfectly aligned spins.
This initial spin locations are given by `positions` (Nx3 matrix or sequence of vectors of size 3).
Alternatively the number of spins can be given in which case the spins are randomly distributed in a 1x1x1 mm box centered on the origin.
"""
struct Snapshot{T<:AbstractFloat}
    spins :: AbstractVector{Spin{T}}
    time :: T
    Snapshot(spins :: AbstractVector{Spin{T}}, time=0.) where {T} = new{T}(spins, time)
end

function Snapshot(positions :: AbstractMatrix{<:Real}; time :: Real=0., kwargs...) 
    @assert size(positions, 2) == 3
    Snapshot(map(p -> Spin(position=p; kwargs...), eachrow(positions)), time)
end
function Snapshot(positions :: AbstractVector{<:AbstractVector{<:Real}}; time :: Real=0., kwargs...) 
    Snapshot(map(p -> Spin(position=p; kwargs...), positions), time)
end
Snapshot(nspins :: Int; kwargs...) = Snapshot(rand(nspins, 3) * 1000 - 500; kwargs...)

"""
    time(snapshot)
    time(sequence_component)
    time(sequence, sequence_index)

Returns the time in milliseconds that a snapshot was taken or that a sequence component will have effect.
"""
Base.time(s :: Snapshot) = s.time

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

"""
Represents the positions and multiple orientations of multiple [`Spin`](@ref) objects at a specific `time`.

This object keeps track of a different spin orientation for each sequence.
The orientation for a single sequence can be extracted using [`get_sequence`](@ref), which will return a [`Snapshot`](@ref).

# Useful constructors
    MultiSnapshot(snap::Snapshot, nsequences::Integer)

This replicates the orientation defined in `snap` to be the starting point for all `n_sequences` sequences.

    MultiSnapshot(positions, nsequences::Integer; time=0., longitudinal=1., transverse=0., phase=0.)
    MultiSnapshot(nspins :: Integer, nsequences::Integer; time=0., longitudinal=1., transverse=0., phase=0.)

These create a new `MultiSnapshot` from scratch in the way described in [`Snapshot`](@ref).
"""
struct MultiSnapshot{N, T<:AbstractFloat}
    # datatype T, N sequences, M spins
    spins :: AbstractVector{MultiSpin{N, T}}
    time :: T
    MultiSnapshot(spins::AbstractVector{MultiSpin{N,T}}, time::Real=0.) where {N,T} = new{N,T}(spins, time)
end
MultiSnapshot(snap :: Snapshot, nsequences::Integer) = MultiSnapshot([MultiSpin(spin, nsequences) for spin in snap.spins], snap.time)
MultiSnapshot(snap, nsequences::Integer; kwargs...) = MultiSnapshot(Snapshot(snap; kwargs...), nsequences)
Base.time(s :: MultiSnapshot) = s.time
get_sequence(snap::MultiSnapshot, index) = Snapshot(get_sequence.(snap.spins, index), snap.time)

for param in (:longitudinal, :transverse, :phase, :vector)
    @eval $param(snap :: MultiSnapshot{N}) where {N} = SVector[N]([get_sequence($param(snap), idx) for idx in 1:N])
end