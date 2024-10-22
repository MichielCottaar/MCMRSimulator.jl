"""
Defines the functions that run the actual simulation:
- [`readout`](@ref): get total signal or [`Snapshot`](@ref) at any [`Readout`](@ref) objects in the sequences.
- [`evolve`](@ref): Return a single [`Snapshot`](@ref) with the state of the simulation at a given time. This snapshot can be used as initialisation for further runs.

All of these functions call [`evolve_to_time`](@ref) under the hood to actually run the simulation.
"""
module Evolve
import StaticArrays: SVector, MVector, StaticVector
import LinearAlgebra: norm, ⋅
import MRIBuilder: Sequence, variables, B0
import MRIBuilder.Components: InstantGradient, InstantPulse
import Rotations
import Bessels: besseli0
import ..SequenceParts: SequencePart, MultSequencePart, InstantSequencePart, iter_parts, nreadouts_per_TR, first_TR_with_all_readouts, IndexedReadout
import ..Methods: get_time
import ..Spins: @spin_rng, Spin, Snapshot, stuck, SpinOrientationSum, get_sequence, orientation, SpinOrientation
import ..Simulations: Simulation, _to_snapshot
import ..Relax: relax!
import ..Properties: GlobalProperties, stick_probability
import ..Subsets: Subset, get_subset
import ..Geometries.Internal: Reflection, detect_intersection, empty_intersection, has_intersection, surface_relaxation, permeability, surface_density, direction, previous_hit, dwell_time, empty_reflection

"""
Supertype for any Readout accumulator.
"""
abstract type ReadoutAccumulator end
abstract type SingleAccumulator <: ReadoutAccumulator end

struct TotalSignalAccumulator <: SingleAccumulator
    nspins :: Ref{Int}
    as_vector :: MVector{3, Float64}
    nwrite :: Ref{Int}
    TotalSignalAccumulator() = new(
        Ref(0),
        zero(MVector{3, Float64}),
        Ref(0.),
        Ref(0)
    )
end

struct SnapshotAccumulator <: SingleAccumulator
    spins :: Vector{Spin{1}}
    time :: Ref{Float64}
    nwrite :: Ref{Int}
    SnapshotAccumulator() = new(
        Spin{1}[],
        Ref(0.),
        Ref(0)
    )
end

struct FillerAccumulator <: SingleAccumulator
    nwrite :: Ref{Int}
end

struct MultipleAccumulator{S, T<:ReadoutAccumulator} <: ReadoutAccumulator
    accumulators :: Dict{Int, T}
    flatten :: bool
end

struct SubsetAccumulator{T<:ReadoutAccumulator}
    geometry :: FixedGeometry
    subsets :: Vector{<:Subset}
    accumulators :: Vector{T}
    flatten :: bool
end

function readout!(acc::TotalSignalAccumulator, spins::Vector{Spin{1}}, ::IndexedReadout)
    acc.nspins[] += length(spins)
    acc.as_vector .+= orientation(spins)
    acc.nwrite[] += 1
end

function readout!(acc::SnapshotAccumulator, spins::Vector{Spin{1}}, index::IndexedReadout)
    append!(acc.spins, deepcopy(spins))
    if iszero(acc.time[])
        acc.time[] = index.time
    else
        @assert acc.time[] == index.time
    end
    acc.nwrite[] += 1
end

function readout!(::FillerAccumulator, spins::Vector{Spin{1}}, ::IndexedReadout)
    error("Should not write to a FillerAccumulator!")
end

function readout!(acc::MultipleAccumulator{S}, spins::Vector{Spin{1}}, index::IndexedReadout) where {S}
    use_index = getproperty(index, S)
    if use_index in keys(acc.accumulators)
        readout!(acc.accumulators[use_index], spins, simulation, index)
    end
end

function readout!(acc::SubsetAccumulator, spins::Vector{Spin{1}}, index::IndexedReadout)
    for (subset, accumulator) in zip(acc.subsets, acc.accumulators)
        readout!(accumulator, get_subset(spins, simulation, subset), acc.geometry, index)
    end
end


"""
    has_nwrite(accumulator, target)

Returns true if the number of times that each of the accumulators has been written matches `target`.
"""
has_nwrite(acc::SingleAccumulator, nwrite::Int) = acc.nwrite[] == nwrite
has_nwrite(acc::Union{SubsetAccumulator, MultipleAccumulator}, nwrite::Int) = all(nspins.(acc.accumulators, nwrite))
has_nwrite(::FillerAccumulator, nwrite::Int) = true

"""
    fix_accumulator(accumulator, spins)

Return the accumulated readout results as an array.
"""
fix_accumulator(acc::TotalSignalAccumulator) = SpinOrientationSum(
    SpinOrientation(acc.as_vector),
    acc.nspins
)
fix_accumulator(acc::TotalSignalAccumulator) = Snapshot(acc.spins, acc.as_vector, acc.time)
fix_accumulator(::FillerAccumulator) = nothing
fix_accumulator(acc::Union{SubsetAccumulator, MultipleAccumulator}) = acc.flatten ? fix_accumulator(acc.accumulators[1]) : stack(fix_accumulator.(acc.accumulators))


function base_accumulator(; return_snapshot=false)
    return return_snapshot ? SnapshotAccumulator() : TotalSignalAccumulator()
end

function subset_accumulator(; subset=Subset(), kwargs...)
    if subset isa Subset
        return SubsetAccumulator([subset], [base_accumulator(; kwargs...)], true)
    else
        return SubsetAccumulator(subset, [base_accumulator(; kwargs...) for _ in 1:length(subset)], false)
    end
end

function TR_accumulator(sequence, current_time, start_TR; skip_TR=nothing, nTR=nothing, kwargs...)
    if isnothing(nTR) && isnothing(skip_TR)
        # sequence will not be repeated
        return subset_accumulator(; kwargs...)
    else
        flatten = isnothing(nTR)
        if flatten
            nTR = 1
        end
        use_start_TR = isnothing(skip_TR) ? start_TR : (start_TR + skip_TR)
        return MultipleAccumulator{:TR}(
            Dict(TR => subset_accumulator(; kwargs...) for TR in use_start_TR:(use_start_TR+nTR-1)),
            flatten
        )
    end
end

function readout_accumulator(sequence, current_time; readout_times=nothing, kwargs...)
    if isnothing(readout_times)
        nreadout = nreadouts_per_TR(sequence)
        flatten = isone(nreadout)
        if iszero(nreadout)
            error("Readout times need to be explicitly set when simulating a sequence without readouts.")
        end
    elseif readout_times isa Number
        flatten = true
        nreadout = 1
    else
        flatten = false
        nreadout = length(readout_times)
        if iszero(nreadout)
            error("Readout times cannot be set to an empty vector.")
        end
    end

    start_TR = first_TR_with_all_readouts(sequence, current_time; readout_times=readout_times)

    if flatten
        return TR_accumulator(sequence, current_time, start_TR; kwargs...)
    else
        return MultipleAccumulator{:readout}(
            Dict(readout => TR_accumulator(sequence, current_time, start_TR; kwargs...) for readout in 1:length(readout_times)),
            false
        )
    end
end

function sequence_accumulator(simulation, current_time; kwargs...)
    flatten = simulation.flatten
    return MultipleAccumulator{:sequence}(
        Dict(index => readout_accumulator(sequence, current_time; kwargs...) for (index, sequence) in enumerate(simulation.sequences)),
        flatten
    )
end








"""
    readout(spins, simulation[, readout_times]; bounding_box=<1x1x1 mm box>, skip_TR=0, nTR=1, return_snapshot=false, subset=<all>)

Evolves a set of spins through the [`Simulation`](@ref).
Returns the total signal or a full [`Snapshot`](@ref) at every readout time in the simulated sequences over one or more repetition times (TRs).

# Positional arguments:
- `spins`: Number of spins to simulate or an already existing [`Snapshot`](@ref).
- `simulation`: [`Simulation`](@ref) object defining the environment, scanner, and sequence(s).
- `times` (optional): time of the readouts relative to the start of the TR (in ms). If not provided, the times of any [`Readout`](@ref) objects in the sequence will be used.

# Keyword arguments:
- `bounding_box`: size of the voxel in which the spins are initiated in um (default is 1000, corresponding to a 1x1x1 mm box centered on zero). Can be set to a [`BoundingBox`](@ref) object for more control.
- `skip_TR`: Number of repetition times to skip before starting the readout. 
    Even if set to zero (the default), the simulator will still skip the current TR before starting the readout 
    if the starting snapshot is from a time past one of the sequence readouts.
- `nTR`: number of TRs for which to store the output
- `return_snapshot`: set to true to output the state of all the spins as a [`Snapshot`](@ref) at each readout instead of a [`SpinOrientationSum`](@ref) with the total signal.
- `subset`: Return the signal/snapshot for a subset of all spins. Can be set to a single or a vector of [`Subset`](@ref) objects. If set to a vector, this will add an attional dimension to the output.

# Returns
The function returns an up to 3-dimensional (KxLxMxN) array, with the following dimensions:
- `K`: the number of sequences. This dimension is not included if the simulation only contains a single sequencen (and this single sequence is not passed into the [`Simulation`](@ref) as a vector).
- `L`: the number of readout times with a single TR. This dimension is skipped if the `readout_times` is set to a scalar number. This dimension might contain `nothing`s for sequences that contain fewer [`Readout`](@ref) objects than the maximum (`M`).
- `M`: the number of TRs (controlled by the `nTR` keyword). If `nTR` is not explicitly set by the user, this dimension is skipped.
- `N`: the number of subsets (controlled by the `subset` keyword). If `subset` is set to a single value (<all> by default), this dimension is skipped.
By default each element of this matrix is either a [`SpinOrientationSum`](@ref) with the total signal.
If `return_snapshot=true` is set, each element is the full [`Snapshot`](@ref) instead.
"""
readout(spins, simulation::Simulation, new_readout_times=nothing; bounding_box=500, kwargs...) = readout(_to_snapshot(spins, simulation, bounding_box), simulation, new_readout_times; kwargs...)

function readout(snapshot::Snapshot{N}, simulation::Simulation{N}, new_readout_times=nothing; noflatten=false, skip_TR=nothing, nTR=nothing, kwargs...) where {N}
    accumulator = sequence_accumulator(simulation, snapshot.time; skip_TR=skip_TR, nTR=nTR, kwargs...)
    repeat = Val(!isnothing(skip_TR) || !isnothing(nTR))
    run_readout!(snapshot, simulation, accumulator, repeat)
    return fix_accumulator(accumulator)
end

"""
    run_readout!(snapshot, simulation, accumulators)

Runs the simulation starting from `snapshot` and filling the `accumulators`.
"""
function run_readout!(snapshot::Snapshot{N}, simulation::Simulation{N}, accumulator::ReadoutAccumulator, repeat::Val, target_nwrite=1)
    for part in iter_parts(simulation.sequences, current_time, repeat, simulation.timestep)
        process_sequence_step!(snapshot.spins, simulation, part, B0s, accumulator)
        if has_nwrite(accumulator, target_nwrite)
            return
        end
    end
    error("Simulation ended before all readout accumulators were filled.")
end

"""
    evolve(snapshot, simulation[, new_time]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) to a new time.
Returns a [`Snapshot`](@ref) at the new time, which can be used as a basis for further simulation.
If undefined `new_time` will be set to the start of the next TR.
"""
function evolve(spins, simulation::Simulation{N}, new_time=nothing; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    if isnothing(new_time)
        if iszero(N)
            error("Simulation time needs to be set explicitly for simulations that do not contain any sequences.")
        end
        first_TR = variables.duration(simulation.sequences[1])
        if !all(variables.duration(s) ≈ first_TR for s in simulation.sequences)
            error("Cannot evolve snapshot for a single TR, because the simulation contains sequences with different TRs. Please set a `new_time` explicitly.")
        end
        new_time = (div(nextfloat(snapshot.time), first_TR, RoundDown) + 1) * first_TR
    end
    evolve_to_time(snapshot, simulation, Float64(new_time))
end

"""
    evolve_to_time(snapshot, simulation, new_time)

Evolves the full [`Snapshot`](@ref) through the [`Simulation`](@ref) to the given `new_time`.
Multi-threading is used to evolve multiple spins in parallel.
This is used internally when calling any of the snapshot evolution methods (e.g., [`evolve`](@ref)).
"""
function evolve_to_time(snapshot::Snapshot{N}, simulation::Simulation{N}, new_time::Float64) where {N}
    current_time::Float64 = snapshot.time
    if new_time < current_time
        error("New requested time ($(new_time)) is less than current time ($(snapshot.time)). Simulator does not work backwards in time.")
    end
    if new_time == current_time
        return snapshot
    end
    B0s = map(B0, simulation.sequences)

    for part in iter_parts(simulation.sequences, current_time, new_time, simulation.timestep)
        process_sequence_step!(snapshot.spins, simulation, part, B0s)
    end
    return Snapshot(snapshot.spins, new_time)
end

process_sequence_step!(spins, simulation, sequence_part::MultSequencePart, B0s, accumulator) = draw_step!(spins, simulation, sequence_part, B0s)
process_sequence_step!(spins, simulation, sequence_part::InstantSequencePart, B0s, accumulator) = apply_instants!(spins, sequence_part, accumulator)

"""
    draw_step!(spin(s), simulation, mult_sequence_part, B0s)

Updates the spin based on a random movement through the given geometry for a given `timestep`:
- draws the next location of the particle after `timestep` with given `simulation.diffusivity`.  
  This displacement will take into account the obstructions in `simulation.geometry`.
- The spin orientation will be affected by relaxation (see [`relax!`](@ref)) and potentially by magnetisation transfer during collisions.
"""
function draw_step!(spins::Vector{Spin{N}}, simulation::Simulation{N}, sequence_part::MultSequencePart{N}, B0s::StaticVector{N, Float64}) where {N}
    Threads.@threads for spin in spins
        draw_step!(spin, simulation, sequence_part, B0s)
    end
end


function draw_step!(spin::Spin{N}, simulation::Simulation{N}, parts::MultSequencePart{N}, B0s::StaticVector{N, Float64}, test_new_pos=nothing) where {N}
    if ~isnothing(test_new_pos)
        all_positions = [spin.position]
    end
    timestep = parts.duration
    if iszero(timestep)
        if ~isnothing(test_new_pos)
            return all_positions
        else
            return
        end
    end
    current_pos = spin.position
    is_stuck = stuck(spin)
    fraction_timestep = 0.

    found_solution = false
    @spin_rng spin begin
        if !stuck(spin)
            if isnothing(test_new_pos)
                rand_base_vec = randn(SVector{3, Float64})
                displacement = rand_base_vec .* sqrt(2 * simulation.diffusivity * timestep)
                new_pos = current_pos + displacement
            else
                new_pos = SVector{3, Float64}(test_new_pos)
            end
            reflection = Reflection(norm(new_pos - current_pos) / sqrt(2 * simulation.diffusivity * timestep))
            phit = (0, 0, false)
        end
        for _ in 1:1000000
            if is_stuck
                td = dwell_time(simulation.geometry, spin.reflection)
                fraction_stuck = -log(rand()) * td / timestep
                relax!(spin, spin.reflection.inside, simulation, parts, fraction_timestep, min(one(Float64), fraction_timestep + fraction_stuck), B0s)
                fraction_timestep += fraction_stuck
                if fraction_timestep >= 1
                    found_solution = true
                    break
                end
                phit = previous_hit(spin.reflection)
                reflection = spin.reflection
                displacement = direction(reflection, (1 - fraction_timestep) * timestep, simulation.diffusivity)
                new_pos = current_pos .+ displacement

                spin.reflection = empty_reflection
                is_stuck = false
            end

            collision = detect_intersection(
                simulation.geometry,
                current_pos,
                new_pos,
                phit,
            )

            use_distance = collision === empty_intersection ? 1. : max(prevfloat(collision.distance), 0.)
            if collision === empty_intersection
                next_fraction_timestep = 1.
                collision_pos = new_pos
            else
                next_fraction_timestep = fraction_timestep + (1 - fraction_timestep) * use_distance
                collision_pos = map((p1, p2) -> use_distance * p1 + (1 - use_distance) * p2, new_pos, current_pos)
            end

            # spin relaxation
            relax!(spin, collision_pos, simulation, parts, fraction_timestep, next_fraction_timestep, B0s)

            if ~has_intersection(collision)
                spin.position = new_pos
                found_solution = true
                break
            end

            relaxation = surface_relaxation(simulation.geometry, collision)
            if ~iszero(relaxation)
                collision_attenuation = exp(-sqrt(timestep) * relaxation)
                for orientation in spin.orientations
                    orientation.transverse *= collision_attenuation
                end
            end

            normed_rate = sqrt(timestep) * permeability(simulation.geometry, collision)
            permeability_prob = isinf(normed_rate) ? 1. : (1. - exp(-normed_rate) * besseli0(normed_rate))
            passes_through = isone(permeability_prob) || !(iszero(permeability_prob) || rand() > permeability_prob)
            reflection = Reflection(collision, new_pos - current_pos, reflection.ratio_displaced, 
                reflection.time_moved + (1 - fraction_timestep) * use_distance * timestep, 
                reflection.distance_moved + norm(new_pos - current_pos) * use_distance, 
                passes_through
            )
            current_pos = spin.position = collision_pos
            if ~isnothing(test_new_pos)
                push!(all_positions, current_pos)
            end

            sd = surface_density(simulation.geometry, collision)
            if !iszero(sd) && rand() < stick_probability(sd, dwell_time(simulation.geometry, collision), simulation.diffusivity, timestep)
                spin.reflection = reflection
                is_stuck = true
            else
                new_pos = current_pos .+ direction(reflection, (1 - next_fraction_timestep) * timestep, simulation.diffusivity)
            end

            phit = previous_hit(reflection)
            fraction_timestep = next_fraction_timestep
        end
    end
    if !found_solution
        error("Bounced single particle for 1000000 times in single step; terminating!")
    end
    if ~isnothing(test_new_pos)
        push!(all_positions, spin.position)
        return all_positions
    end
end

"""
    apply_instants!(spins, instants, accumulator)

Apply a set of `N` instants to a vector of spins for each of the `N` sequences being simulated.

Each instant can be:
- `nothing`: do nothing
- `MRIBuilder.Components.InstantPulse`: apply RF pulse rotation
- `MRIBuilder.Components.InstantGradient`: add phase corresponding to gradient
- `IndexedReadout`: add the current snapshot to `accumulator`
"""
apply_instants!(spins::Vector{Spin{N}}, instants::InstantSequencePart{N}, accumulator) where {N} = apply_instants!(spins, instants.instants, accumulator)
function apply_instants!(spins::Vector{Spin{N}}, instants::StaticVector{N}, accumulator) where {N}
    for i in 1:N
        apply_instants!(spins, i, instants[i], accumulator)
    end
end

apply_instants!(spins::Vector{Spin{N}}, instants::StaticVector{N, Nothing}, _) where {N} = nothing
apply_instants!(spins::Vector{<:Spin}, index::Int, ::Nothing, _) = nothing

function apply_instants!(spins::Vector{<:Spin}, index::Int, grad::InstantGradient, _)
    Threads.@threads for spin in spins
        new_phase = rad2deg(spin.position ⋅ variables.qvec(grad))
        spin.orientations[index].phase += new_phase
    end
end

function apply_instants!(spins::Vector{<:Spin{1}}, index::Int, readout::IndexedReadout, accumulator::ReadoutAccumulator)
    readout!(accumulator, spins, readout)
end

function apply_instants!(spins::Vector{<:Spin}, index::Int, pulse::InstantPulse, _)
    rotation = Rotations.RotationVec(
        deg2rad(pulse.flip_angle) * cosd(pulse.phase),
        deg2rad(pulse.flip_angle) * sind(pulse.phase),
        0.
    )
    Threads.@threads for spin in spins
        orient = spin.orientations[index]
        new_orient = SpinOrientation(rotation * orientation(orient))
        orient.longitudinal = new_orient.longitudinal
        orient.transverse = new_orient.transverse
        orient.phase = new_orient.phase
    end
end

end