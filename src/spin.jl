"""
    norm_angle(angle)

Normalises an angle in degrees, so that it is between it is in the range (-180, 180]
"""
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
mutable struct SpinOrientation
    longitudinal :: Float
    transverse :: Float
    phase :: Float
    SpinOrientation(longitudinal, transverse, phase) = new(Float(longitudinal), Float(transverse), Float(phase))
end

SpinOrientation(orientation::AbstractVector) = SpinOrientation(PosVector(orientation))
SpinOrientation(vector::PosVector) = SpinOrientation(
    vector[3],
    sqrt(vector[1] * vector[1] + vector[2] * vector[2]),
    rad2deg(atan(vector[2], vector[1])),
)

Base.show(io::IO, orient::SpinOrientation) = print(io, "SpinOrientation(longitudinal=$(longitudinal(orient)), transverse=$(transverse(orient)), phase=$(phase(orient))°)")

"""
Spin particle with a position and `nsequences` spin orientations (stored as [`SpinOrientation`](@ref)).

A random number generator is stored in the `Spin` object as well, which will be used for evolving the spin into the future in a reproducible manner.

# Constructors
    Spin(;nsequences=1, position=[0, 0, 0], longitudinal=1., transverse=0., phase=0.)

Creates a new spin with `nsequences` identical spin orientations (given by `longitudinal`, `transverse`, and `phase` flags).
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
mutable struct Spin{N}
    position :: PosVector
    orientations :: SVector{N, SpinOrientation}
    stuck_to :: Reflection
    rng :: FixedXoshiro
    function Spin(position::AbstractArray{<:Real}, orientations::AbstractArray{SpinOrientation}, stuck_to=empty_reflection, rng::FixedXoshiro=FixedXoshiro()) 
        new{length(orientations)}(PosVector(position), SVector{length(orientations)}(deepcopy.(orientations)), stuck_to, rng)
    end
end

function Spin(;nsequences=1, position=zero(SVector{3,Float}), longitudinal=1., transverse=0., phase=0., stuck_to=empty_reflection, rng=FixedXoshiro()) 
    base = Spin(SVector{3, Float}(position), SVector{1}(SpinOrientation(longitudinal, transverse, phase)), stuck_to, rng)
    return nsequences == 1 ? base : Spin(base, nsequences)
end
Spin(reference_spin::Spin{1}, nsequences::Int) = Spin(reference_spin.position, repeat(reference_spin.orientations, nsequences), reference_spin.stuck_to, reference_spin.rng)

show_helper(io::IO, spin::Spin{0}) = print(io, "with no magnetisation information)")
show_helper(io::IO, spin::Spin{1}) = print(io, "with $(repr(spin.orientations[1], context=io)))")
show_helper(io::IO, spin::Spin{N}) where {N} = print(io, "with magnetisations for $N sequences)")
function Base.show(io::IO, spin::Spin) 
    if stuck(spin)
        print(io, "stuck ")
    end
    print(io, "Spin(position=$(spin.position) ")
    show_helper(io, spin)
end
stuck(spin::Spin) = Reflection(spin) !== empty_reflection
Reflection(spin) = spin.stuck_to
Collision(spin) = Reflection(spin).collision
stuck_to(spin) = Collision(spin).properties

macro spin_rng(spin, expr)
    return quote
        local old_rng_state = copy(Random.TaskLocalRNG())
        copy!(Random.TaskLocalRNG(), $(esc(spin)).rng)
        result = $(esc(expr))
        $(esc(spin)).rng = FixedXoshiro(copy(Random.TaskLocalRNG()))
        copy!(Random.TaskLocalRNG(), old_rng_state)
        result
    end
end

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
phase(o :: SpinOrientation) = norm_angle(o.phase)
orientation(o :: SpinOrientation) = PosVector(
    o.transverse * cosd(o.phase),
    o.transverse * sind(o.phase),
    o.longitudinal
)

for param in (:longitudinal, :transverse, :phase, :orientation)
    @eval $param(s :: Spin{1}) = $param(s.orientations[1])
    @eval $param(s :: Spin) = map($param, s.orientations)
end

"""
    get_sequence(spin, sequence_index)
    get_sequence(snapshot, sequence_index)

Extracts the spin orientation corresponding to a specific sequence, where the sequence index uses the order in which the sequences where provided in the [`Simulation`](@ref).
"""
get_sequence(spin::Spin, index) = Spin(spin.position, SVector{1}([spin.orientations[index]]), spin.stuck_to, spin.rng)


"""
    position(s::Spin)

Returns the position of the spin particle as a vector of length 3.
"""
position(s :: Spin) = s.position

isinside(bb::BoundingBox, spin::Spin) = isinside(bb, position(spin))
isinside(other, spin::Spin) = isinside(other, position(spin), Collision(spin))
inside_MRI_properties(geom::Geometry, spin::Spin, global_props::MRIProperties) = inside_MRI_properties(geom, position(spin), global_props)

"""
Represents the positions and orientations of multiple [`Spin`](@ref) objects at a specific `time`.

Note that times are in milliseconds and positions in micrometer. 
The equilibrium longitudinal spin (after T1 relaxation) is always 1.

# Useful constructors
    Snapshot(positions; time=0., longitudinal=1., transverse=0., phase=0., nsequences=1)
    Snapshot(nspins[, bounding_box[, geometry, default_surface_density]]; time=0., longitudinal=1., transverse=0., phase=0., nsequences=1)
    Snapshot(nspins, simulation[, bounding_box; time=0., longitudinal=1., transverse=0., phase=0., nsequences=1)

Creates a new Snapshot at the given `time` with spins initialised for simulating `nsequences` sequences.
All spins will start out in equilibrium, but that can be changed using the `longitudinal`, `transverse`, and/or `phase` flags.
This initial spin locations are given by `positions` (Nx3 matrix or sequence of vectors of size 3).
Alternatively the number of spins can be given in which case the spins are randomly distributed in the given `bounding_box` (default: 1x1x1 mm box centered on origin).
The bounding_box can be a [`BoundingBox`](@ref) object, a tuple with the lower and upper bounds (i.e., two vectors of length 3) or a number `r` (resulting in spins filling a cube from `-r` to `+r`)

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
struct Snapshot{N} <: AbstractVector{Spin{N}}
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

function Snapshot(nspins::Integer, bounding_box=500, geometry=nothing, default_surface_density=zero(Float); time::Real=0., kwargs...)
    bounding_box = BoundingBox(bounding_box)
    sz = (upper(bounding_box) - lower(bounding_box))
    free_spins = map(i->Spin(; position=rand(PosVector) .* sz .+ lower(bounding_box), kwargs...), 1:nspins)
    if isnothing(geometry)
        return Snapshot(free_spins, time)
    end
    geometry = Geometry(geometry)
    volume = prod(sz)
    density = nspins / volume
    stuck_spins = random_surface_spins(geometry, bounding_box, density, default_surface_density; kwargs...)
    if length(stuck_spins) == 0
        spins = free_spins
    else
        spins = Random.shuffle(vcat(free_spins, stuck_spins))[1:nspins]
    end
    return Snapshot(spins, time)
end

Snapshot(nspins :: Int; kwargs...) = Snapshot(rand(nspins, 3) .* 1000 .- 500; kwargs...)

Base.show(io::IO, snap::Snapshot{1}) = print(io, "Snapshot($(length(snap)) spins with total magnetisation of $(repr(SpinOrientation(snap), context=io)) at t=$(get_time(snap))ms)")
Base.show(io::IO, snap::Snapshot{N}) where {N} = print(io, "Snapshot($(length(snap)) spins with magnetisations for $N sequences at t=$(get_time(snap))ms)")

"""
    get_time(snapshot)
    get_time(sequence_component)
    get_time(sequence, sequence_index)

Returns the time in milliseconds that a snapshot was taken or that a sequence component will have effect.
"""
get_time(s :: Snapshot) = s.time

function orientation(s :: Snapshot)
    sum(orientation, s.spins)
end
SpinOrientation(s :: Snapshot) = SpinOrientation(orientation(s))

for param in (:longitudinal, :transverse, :phase)
    @eval $param(s :: Snapshot) = $param(SpinOrientation(s))
end

# Abstract Vector interface (following https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)
Base.size(s::Snapshot) = size(s.spins)
Base.getindex(s::Snapshot, index::Int) = s.spins[index]
Base.getindex(s::Snapshot, index) = Snapshot(s.spins[index], s.time)
Base.length(s::Snapshot) = length(s.spins)
Base.iterate(s::Snapshot) = iterate(s.spins)
Base.iterate(s::Snapshot, state) = iterate(s.spins, state)
Base.IndexStyle(::Type{Snapshot}) = IndexLinear()

"""
    position.(s::Snapshot)

Returns all the positions of the spin particles as a vector of length-3 vectors.
"""
position(s::Snapshot) = position.(s.spins)

Snapshot(snap :: Snapshot{1}, nsequences::Integer) = Snapshot([Spin(spin, nsequences) for spin in snap.spins], snap.time)
get_sequence(snap::Snapshot, index) = Snapshot(get_sequence.(snap.spins, index), snap.time)
isinside(something, snapshot::Snapshot) = [isinside(something, spin) for spin in snapshot]


project(something, v::AbstractVector) = project(something, SVector{length(v)}(v))
project(something, spin::Spin) = project(something, position(spin))

function project(something, snap::Snapshot)
    Snapshot(
        [Spin(project(something, s), s.orientations) for s in snap.spins],
        snap.time
    )
end
