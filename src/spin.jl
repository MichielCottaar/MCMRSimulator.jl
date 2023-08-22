"""
Types:
- [`Snapshot`](@ref)
- [`Spin`](@ref)
- [`SpinOrientation`](@ref)
- [`FixedXoshiro`](@ref)

Methods:
- [`longitudinal`](@ref)
- [`transverse`](@ref)
- [`phase`](@ref)
- [`orientation`](@ref)
- [`isinside`](@ref)
- [`update_isinside`](@ref)
- [`stuck`](@ref)
- [`stuck_to`](@ref)
- [`get_sequence`](@ref)
"""
module Spins

import Random
import StaticArrays: SVector
import LinearAlgebra: ⋅, norm
import ..Geometries.Internal: 
    Reflection, empty_reflection, has_intersection,
    BoundingBox, lower, upper, random_surface_positions,
    FixedGeometry, FixedObstruction, FixedObstructionGroup,
    isinside, R1, R2, off_resonance, has_hit, previous_hit
import ..Methods: get_time, norm_angle, phase, project
import ..Properties: GlobalProperties
import ..Geometries: fix

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
    longitudinal :: Float64
    transverse :: Float64
    phase :: Float64
    SpinOrientation(longitudinal, transverse, phase) = new(Float64(longitudinal), Float64(transverse), Float64(phase))
end

SpinOrientation(orientation::AbstractVector) = SpinOrientation(SVector{3, Float64}(orientation))
SpinOrientation(vector::SVector{3, Float64}) = SpinOrientation(
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
    position :: SVector{3, Float64}
    orientations :: SVector{N, SpinOrientation}
    reflection :: Reflection
    inside :: Vector{Tuple{Int, Int}}
    rng :: FixedXoshiro
    function Spin(position::AbstractArray{<:Real}, orientations::AbstractArray{SpinOrientation}, reflection=empty_reflection, rng::FixedXoshiro=FixedXoshiro()) 
        new{length(orientations)}(SVector{3, Float64}(position), SVector{length(orientations)}(deepcopy.(orientations)), reflection, Tuple{Int, Int}[], rng)
    end
end

function Spin(;nsequences=1, position=zero(SVector{3,Float64}), longitudinal=1., transverse=0., phase=0., reflection=empty_reflection, rng=FixedXoshiro()) 
    base = Spin(SVector{3, Float64}(position), SVector{1}(SpinOrientation(longitudinal, transverse, phase)), reflection, rng)
    return nsequences == 1 ? base : Spin(base, nsequences)
end
Spin(reference_spin::Spin{1}, nsequences::Int) = Spin(reference_spin.position, repeat(reference_spin.orientations, nsequences), reference_spin.reflection, reference_spin.rng)

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

"""
    stuck(spin)

Returns true if the spin is stuck on the surface.
This can be used to filter a [`Snapshot`](@ref) using:
```julia
only_stuck = filter(stuck, snapshot)
only_free = filter(s -> !stuck(s), snapshot)
```
"""
stuck(spin::Spin) = has_intersection(spin.reflection)

"""
    stuck_to(spin)

Return the indices of the obstruction the spin is stuck to.
Will return (0, 0) for a free spin.
"""
function stuck_to(spin::Spin)
    return has_hit(spin.reflection)
end

"""
    stuck_to(spin, geometry)

Return the internal representation of the obstruction the spin is stuck to.
Raises an error if the spin is free.
"""
function stuck_to(spin::Spin, geometry)
    (outer_index, inner_index) = has_hit(spin)
    if iszero(outer_index)
        error("Free spin is not stuck to any obstructions.")
    end
    return geometry[outer_index][inner_index]
end


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
    longitudinal(spin)
    longitudinal(snapshot)

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
orientation(o :: SpinOrientation) = SVector{3, Float64}(
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
get_sequence(spin::Spin, index) = Spin(spin.position, SVector{1}([spin.orientations[index]]), spin.reflection, spin.rng)


"""
    position(s::Spin)

Returns the position of the spin particle as a vector of length 3.
"""
position(s :: Spin) = s.position

"""
    isinside([geometry, ]spin)

Returns vector of obstructions that the spin is inside.
If `geometry` is not provided, will return a vector of indices instead.
If a non-fixed `geometry` is provided, will return the number of obstructions that the spin is inside.
"""
isinside(geometry, position::AbstractVector) = isinside(geometry, Spin(nsequences=0, position=position))
isinside(geometry, spin::Spin) = length(isinside(fix(geometry), spin))

isinside(geometry::FixedGeometry, spin::Spin) = isinside(geometry, spin.position, spin.reflection)

"""
    property_values(spin, T, symbol, geometry[, global_properties])

Get sequence of values for property `symbol` of type `T` for the given [`Spin`](@ref).
In order this will return:
1. The surface MRI properties (only if the spin is stuck).
2. Any obstructions in the geometry that the spin is inside. These will be applied in the order that they appear in the [`FixedGeometry`](@ref).
3. The global property
Any values of `nothing` will be filtered out.
"""
function property_values(spin::Spin, ::Type{T}, symbol::Symbol, geometry::FixedGeometry, global_properties::GlobalProperties=GlobalProperties()) where {T}
    result = property_values(T, symbol, geometry, spin.inside, stuck_to(spin))
    if hasproperty(global_properties, symbol) && ~isnothing(getproperty(global_properties, symbol))
        push!(result, getproperty(global_properties, symbol))
    end
    return result
end


for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        """
            $($symbol)(spin, geometry, global_properties)
            $($symbol)(position, geometry, global_properties[, stuck_to])
        
        Returns the $($symbol) experienced by the [`Spin`](@ref) given the surface and volume properties of the [`FixedGeometry`](@ref).
        Alternatively, the `position` of the spin can be provided. In this case the [`Reflection`](@ref) should also be returned for a bound spin.
        """
        function $symbol(spin::Spin, geometry::FixedGeometry, global_properties::GlobalProperties=GlobalProperties())
            $symbol(spin, geometry, global_properties, spin.reflection)
        end
        function $symbol(spin_or_pos, geometry, global_properties::GlobalProperties=GlobalProperties())
            $symbol(spin_or_pos, fix(geometry), global_properties)
        end
    end
end

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
    time :: Float64
    Snapshot(spins :: AbstractVector{Spin{N}}, time=0.) where {N} = new{N}(spins, Float64(time))
end

function Snapshot(positions :: AbstractMatrix{<:Real}; time :: Real=0., kwargs...) 
    @assert size(positions, 2) == 3
    Snapshot(map(p -> Spin(; position=p, kwargs...), eachrow(positions)), time)
end
function Snapshot(positions :: AbstractVector{<:AbstractVector{<:Real}}; time :: Real=0., kwargs...) 
    Snapshot(map(p -> Spin(; position=p, kwargs...), positions), time)
end

function Snapshot(nspins::Integer, bounding_box=500, geometry=nothing, default_surface_density=zero(Float64); time::Real=0., kwargs...)
    bounding_box = BoundingBox{3}(bounding_box)
    sz = (upper(bounding_box) - lower(bounding_box))
    free_spins = map(i->Spin(; position=rand(SVector{3, Float64}) .* sz .+ lower(bounding_box), kwargs...), 1:nspins)
    geometry = geometry isa FixedGeometry ? geometry : fix(geometry)
    if iszero(length(geometry))
        return Snapshot(free_spins, time)
    end
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


function random_surface_spins(geometry::FixedGeometry, bounding_box::BoundingBox{3}, volume_density::Number, default_surface_density::Number; nsequences=1, kwargs...)
    spins = Spin{nsequences}[]
    for (position, normal, geometry_index, obstruction_index) in random_surface_positions(geometry, bounding_box, volume_density, default_surface_density)
        inside = Random.rand() > 0.5
        use_normal = inside ? normal : -normal
        direction = Random.randn(SVector{3, Float64})
        if direction ⋅ use_normal < 0
            direction = - direction
        end
        displacement = norm(direction)
        reflection = Reflection(geometry_index, obstruction_index, inside, direction ./ displacement, displacement, 0., 0.)
        push!(spins, Spin(; position=position, reflection=reflection, nsequences=nsequences, kwargs...))
    end
    return spins
end

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

end