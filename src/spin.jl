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

function FixedXoshiro(seed=nothing) 
    if isnothing(seed)
        seed = rand(typemin(UInt64):typemax(UInt64))
    end
    FixedXoshiro(Random.Xoshiro(seed))
end
Random.Xoshiro(rng::FixedXoshiro) = Random.Xoshiro(rng.s0, rng.s2, rng.s2, rng.s3)
function Base.copy!(dst::Random.TaskLocalRNG, src::FixedXoshiro)
    copy!(dst, Random.Xoshiro(src))
    dst
end

"""
    SpinOrientation(longitudinal, transverse, phase)

The spin orientation. Usually created as part of a [`Spin`](@ref) object.

    SpinOrientation(snapshot::Snapshot)

Returns the average spin orientations of all [`Spin`](@ref) objects in the [`Snapshot`](@ref).

This information can be extracted using:
- [`longitudinal`](@ref) to get the spin in the z-direction (equilibrium of 1)
- [`transverse`](@ref) to get the spin in the x-y-plane
- [`phase`](@ref) to get the spin angle in x-y plane (in degrees)
- [`orientation`](@ref) to get the spin orientation as a length-3 vector
"""
struct SpinOrientation
    longitudinal :: Float
    transverse :: Float
    phase :: Float
    SpinOrientation(longitudinal, transverse, phase) = new(Float(longitudinal), Float(transverse), Float(phase))
end

SpinOrientation(orientation::AbstractVector) = SpinOrientation(PosVector(orientation))
SpinOrientation(vector::PosVector) = SpinOrientation(
    vector[3],
    sqrt(vector[1] * vector[1] + vector[2] * vector[2]),
    atan(vector[2], vector[1]),
)

"""
Spin particle with a position and `nsequences` spin orientations (stored as [`SpinOrientation`](@ref)).

A random number generator is stored in the `Spin` object as well, which will be used for evolving the spin into the future in a reproducible manner.

# Constructors
    Spin(;nsequences=1, position=[0, 0, 0], longitudinal=1., transverse=0., phase=0.)

Creates a new spin with `nsequences` identical spin orienations (given by `longitudinal`, `transverse`, and `phase` flags).
The spin will start at given position.

    Spin(reference_spin::Spin{1}, nsequences)

Create a new spin with the same position as `reference_spin` with the orientation of `reference_spin` replicated `nsequences` times.

# Extracting spin information
- [`longitudinal`](@ref) to get the `nsequences` spin magnitudes in the z-direction (equilibrium of 1)
- [`transverse`](@ref) to get the `nsequences` spin magnitudes in the x-y-plane
- [`phase`](@ref) to get the `nsequences` spin angles in x-y plane (in degrees)
- [`orientation`](@ref) to get a (`nsequences`x3) matrix with the spin orientations in 3D space
- [`position`](@ref) to get a length-3 vector with spin location
"""
struct Spin{N}
    position :: PosVector
    orientations :: SVector{N, SpinOrientation}
    rng :: FixedXoshiro
end
function Spin(position::AbstractArray{<:Real}, orientations::AbstractArray{SpinOrientation}, rng::FixedXoshiro=FixedXoshiro()) 
    Spin(PosVector(position), SVector{length(orientations)}(orientations), rng)
end

Spin(;nsequences=1, position=zero(SVector{3,Float}), longitudinal=1., transverse=0., phase=0., rng=FixedXoshiro()) = Spin{nsequences}(SVector{3, Float}(position), SVector{1}(SpinOrientation(longitudinal, transverse, deg2rad(phase))), rng)
Spin(reference_spin::Spin{1}, nsequences::Int) = Spin(reference_spin.position, repeat(reference_spin.orientations, nsequences), reference_spin.rng)

"""
    transverse(spin)
    transverse(snapshot)

Returns the longitudinal magnitude of the spin (i.e., magnitude aligned with the magnetic field) for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available an array of longitudinal values is returned with a value for each sequence.
"""
function longitudinal end

"""
    transverse(spin)
    transverse(snapshot)

Returns the transverse spin (i.e., magnitude in the plane perpendicular to the magnetic field) for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available an array of transverse values is returned with a value for each sequence.
"""
function transverse end

"""
    phase(spin)
    phase(snapshot)

Returns the phase in the x-y plane of the spin for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available  an array of phase values is returned with a value for each sequence.
"""
function phase end

"""
    orientation(spin)
    orientation(snapshot)

Returns the spin orientation as a length-3 vector for a single particle ([`Spin`](@ref)) or averaged across a group of particles in a [`Snapshot`].
When orientations for multiple sequences are available an array of vectors is returned with a value for each sequence.
"""
function orientation end

for param in (:longitudinal, :transverse)
    @eval $param(o :: SpinOrientation) = o.$param
end
phase(o :: SpinOrientation) = norm_angle(rad2deg(o.phase))
orientation(o :: SpinOrientation) = SA[
    o.transverse * cos(o.phase),
    o.transverse * sin(o.phase),
    o.longitudinal
]

for param in (:longitudinal, :transverse, :phase, :orientation)
    @eval $param(s :: Spin{1}) = $param(s.orientations[1])
    @eval $param(s :: Spin) = map($param, s.orientations)
end

"""
    get_sequence(spin, sequence_index)
    get_sequence(snapshot, sequence_index)

Extracts the spin orientation corresponding to a specific sequence, where the sequence index uses the order in which the sequences where provided in the [`Simulation`](@ref).
"""
get_sequence(spin::Spin, index) = Spin(spin.position, SVector{1}([spin.orientations[index]]), spin.rng)


"""
    position(s::Spin)

Returns the position of the spin particle as a vector of length 3.
"""
position(s :: Spin) = s.position

"""
Represents the positions and orientations of multiple [`Spin`](@ref) objects at a specific `time`.

Note that times are in milliseconds and positions in micrometer. 
The equilibrium longitudinal spin (after T1 relaxation) is always 1.

# Useful constructors
    Snapshot(positions; time=0., longitudinal=1., transverse=0., phase=0., nsequences=1)
    Snapshot(nspins::Integer; time=0., longitudinal=1., transverse=0., phase=0., nsequences=1)

Creates a new Snapshot at the given `time` with perfectly longitudinal spins initialised for simulating `nsequences` sequences.
This initial spin locations are given by `positions` (Nx3 matrix or sequence of vectors of size 3).
Alternatively the number of spins can be given in which case the spins are randomly distributed in a 1x1x1 mm box centered on the origin.

    Snapshot(snap::Snapshot{1}, nsequences)

Replicates the positions and orientations for a single sequence in the input snapshot across `nsequences`.

# Extracting summary information
- [`longitudinal`](@ref)(snapshot) to get the `nsequences` spin magnitudes in the z-direction (equilibrium of 1) averaged over all spins
- [`transverse`](@ref)(snapshot) to get the `nsequences` spin magnitudes in the x-y-plane averaged over all spins
- [`phase`](@ref)(snapshot) to get the `nsequences` spin angles in x-y plane (in degrees) averaged over all spins
- [`orientation`](@ref)(snapshot) to get a (`nsequences`x3) matrix with the spin orientations in 3D space averaged over all spins
- [`SpinOrientation`](@ref)(snapshot) to get a `nsequences` vector of [`SpinOrientation`] objects with the average spin orientation across all spins
- [`position`](@ref).(snapshot) to get a the position for each spin in a vector (no averaging applied)

Information for a single sequence can be extracted by calling [`get_sequence`](@ref) first.
"""
struct Snapshot{N}
    spins :: AbstractVector{Spin{N}}
    time :: Float
    Snapshot(spins :: AbstractVector{Spin{N}}, time=0.) where {N} = new{N}(spins, Float(time))
end

function Snapshot(positions :: AbstractMatrix{<:Real}; time :: Real=0., kwargs...) 
    @assert size(positions, 2) == 3
    Snapshot(map(p -> Spin(; position=p, kwargs...), eachrow(positions)), time)
end
function Snapshot(positions :: AbstractVector{<:AbstractVector{<:Real}}; time :: Real=0., kwargs...) 
    Snapshot(map(p -> Spin(; position=p, kwargs...), positions), time)
end
Snapshot(nspins :: Int; kwargs...) = Snapshot(rand(nspins, 3) .* 1000 .- 500; kwargs...)

"""
    time(snapshot)
    time(sequence_component)
    time(sequence, sequence_index)

Returns the time in milliseconds that a snapshot was taken or that a sequence component will have effect.
"""
Base.time(s :: Snapshot) = s.time

function orientation(s :: Snapshot)
    sum(orientation, s.spins)
end
SpinOrientation(s :: Snapshot) = SpinOrientation(orientation(s))

for param in (:longitudinal, :transverse, :phase)
    @eval $param(s :: Snapshot) = $param(SpinOrientation(s))
end
Base.getindex(s::Snapshot, index::Int) = s.spins[index]
Base.getindex(s::Snapshot, index) = Snapshot(s.spins[index], s.time)
Base.length(s::Snapshot) = length(s.spins)
Base.iterate(s::Snapshot) = iterate(s.spins)
Base.iterate(s::Snapshot, state) = iterate(s.spins, state)

"""
    position.(s::Snapshot)

Returns all the positions of the spin particles as a vector of length-3 vectors.
"""
position(s::Snapshot) = position.(s.spins)

Snapshot(snap :: Snapshot{1}, nsequences::Integer) = Snapshot([Spin(spin, nsequences) for spin in snap.spins], snap.time)
get_sequence(snap::Snapshot, index) = Snapshot(get_sequence.(snap.spins, index), snap.time)