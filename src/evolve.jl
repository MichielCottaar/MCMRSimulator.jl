"""
Defines the functions that run the actual simulation:
- [`readout`](@ref): get total signal or [`Snapshot`](@ref) at any [`Readout`](@ref) objects in the sequences.
- [`custom_readout`](@ref): get total signal or [`Snapshot`](@ref) at user-defined times.
- [`evolve`](@ref): Return a single [`Snapshot`](@ref) with the state of the simulation at a given time. This snapshot can be used as initialisation for further runs.

All of these functions call [`evolve_to_time`](@ref) under the hood to actually run the simulation.
"""
module Evolve
import StaticArrays: SVector, MVector
import LinearAlgebra: norm, ⋅
import ..Methods: get_time
import ..Spins: @spin_rng, Spin, Snapshot, stuck, SpinOrientation, get_sequence
import ..Sequences: SequencePart, next_instant, InstantComponent, apply!
import ..Geometries.Internal: 
    FixedGeometry, empty_reflection, detect_intersection, Intersection,
    empty_intersection, direction, Reflection, has_intersection, previous_hit,
    surface_relaxivity, surface_density, dwell_time, permeability, FixedSusceptibility
import ..Simulations: Simulation, _to_snapshot
import ..Relax: relax!, transfer!
import ..Timestep: propose_times
import ..Properties: GlobalProperties, correct_for_timestep, stick_probability

"""
    readout(snapshot, simulation; bounding_box=<1x1x1 mm box>, skip_TR=0, nTR=1)

Evolves the spins in the [`Snapshot`](@ref) through the [`Simulation`](@ref).
Returns the [`Snapshot`](@ref) at every [`Readout`](@ref) in the simulated sequences over one or more repetition times.

The return object depends on whether the simulation was created with a single sequence object or with a vector of sequences.
- For a single sequence object a vector of [`Snapshot`](@ref) objects is returned with a single snapshot for each [`Readout`](@ref) in the sequence during the TRs.
- For a vector of sequences a vector of vectors of [`Snapshot`](@ref) objects is returned. Each element in the outer vector contains the result for a single sequence.


"""
function readout(spins, simulation::Simulation{N}; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)

    readout_times = Float64[]
    per_seq_times = Vector{Float64}[]

    for seq in simulation.sequences
        if length(seq.readout_times) == 0
            seq_times = Float64[]
        else
            current_TR = Int(div(get_time(snapshot), seq.TR, RoundDown))
            time_in_TR = get_time(snapshot) - current_TR * seq.TR

            if minimum(seq.readout_times) < time_in_TR
                current_TR += 1
            end
            seq_times = seq.readout_times .+ (current_TR * seq.TR)
        end
        append!(readout_times, seq_times)
        push!(per_seq_times, seq_times)
    end

    final_snapshots = SVector{N}([Snapshot{1}[] for _ in 1:N])
    for time in sort(unique(readout_times))
        snapshot = evolve_to_time(snapshot, simulation, time)
        for index in 1:N
            if time in per_seq_times[index]
                push!(final_snapshots[index], get_sequence(snapshot, index))
            end
        end
    end
    if simulation.flatten
        return final_snapshots[1]
    else
        return final_snapshots
    end
end

"""
    trajectory(snapshot, simulation, times=[TR]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs at the requested times.
Returns a vector of [`Snapshot`](@ref) objects with the current state of each time in times.
When you are only interested in the signal at each timepoint, use [`signal`](@ref) instead.
"""
function trajectory(spins, simulation::Simulation{N}, times; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    result = Array{typeof(snapshot)}(undef, size(times))
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float64(times[index]))
        result[index] = snapshot
    end
    result
end

"""
    signal(snapshot, simulation, times=[TR]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs the total signal at the requested times.
To get the full snapshot at each timepoint use [`trajectory`](@ref).

The return object depends on whether the simulation was created with a single sequence object or with a vector of sequences.
- For a single sequence object a vector of [`SpinOrientation`](@ref) object with the total signal at each time is returned.
- For a vector of sequences the return type is a (Nt, Ns) matrix of [`SpinOrientation`](@ref) objects for `Nt` timepoints and `Ns` sequences.
"""
function signal(spins, simulation::Simulation{N}, times; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    if simulation.flatten
        result = Array{SpinOrientation}(undef, size(times))
    else
        result = Array{SpinOrientation}(undef, (length(times), N))
    end
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float64(times[index]))
        if simulation.flatten
            result[index] = SpinOrientation(snapshot)
        else
            for seq in 1:N
                result[index, seq] = SpinOrientation(get_sequence(snapshot, seq))
            end
        end
    end
    result
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
        TR = simulation.sequences[1].TR
        if !all(s.TR == TR for s in simulation.sequences)
            error("Cannot evolve snapshot for a single TR, because the simulation contains sequences with different TRs. Please set a `new_time` explicitly.")
        end
        new_time = (div(nextfloat(snapshot.time), TR, RoundDown) + 1) * TR
    end
    evolve_to_time(snapshot, simulation, Float64(new_time))
end

"""
    draw_step!(spin, simulation, sequence_parts, timestep)

Updates the spin based on a random movement through the given geometry for a given `timestep`:
- draws the next location of the particle after `timestep` with given `simulation.diffusivity`.  
  This displacement will take into account the obstructions in `simulation.geometry`.
- The spin orientation will be affected by relaxation (see [`relax!`](@ref)) and potentially by magnetisation transfer during collisions.
"""
function draw_step!(spin :: Spin{N}, simulation::Simulation{N, 0}, parts::SVector{N, SequencePart}, timestep :: Float64) where {N}
    if iszero(timestep)
        return
    end
    relax!(spin, parts, simulation, 0, 1//2)
    @spin_rng spin begin
        spin.position += randn(SVector{3, Float64}) .* sqrt(2 * simulation.diffusivity * timestep)
    end
    relax!(spin, parts, simulation, 1//2, 1)
end

function draw_step!(spin :: Spin{N}, simulation::Simulation{N}, parts::SVector{N, SequencePart}, timestep :: Float64, test_new_pos=nothing) where {N}
    if ~isnothing(test_new_pos)
        all_positions = [spin.position]
    end
    if iszero(timestep)
        if ~isnothing(test_new_pos)
            return all_positions
        else
            return
        end
    end
    current_pos = spin.position
    is_stuck = stuck(spin)
    fraction_timestep = zero(Float64)

    found_solution = false
    @spin_rng spin begin
        if !stuck(spin)
            if isnothing(test_new_pos)
                rand_base_vec = randn(SVector{3, Float64})
                displacement = rand_base_vec .* sqrt(2 * simulation.diffusivity * timestep)
                new_pos = current_pos + displacement
            else
                new_pos = SVector{3}(test_new_pos)
            end
            reflection = Reflection(norm(new_pos - current_pos) / sqrt(2 * simulation.diffusivity * timestep))
            phit = (0, 0, false)
        end
        for _ in 1:10000
            if is_stuck
                td = dwell_time(simulation.geometry, simulation.properties, spin.reflection)
                fraction_stuck = -log(rand()) * td / timestep
                relax!(spin, parts, simulation, fraction_timestep, min(one(Float64), fraction_timestep + fraction_stuck))
                fraction_timestep += fraction_stuck
                if fraction_timestep >= 1
                    found_solution = true
                    break
                end
                phit = previous_hit(spin.reflection)
                reflection = spin.reflection
                displacement = direction(reflection, (1 - fraction_timestep) * timestep, simulation.diffusivity)
                new_pos = current_pos + displacement

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
            next_fraction_timestep = fraction_timestep + (1 - fraction_timestep) * use_distance

            # spin relaxation
            relax_pos_dist = rand() * use_distance
            spin.position = relax_pos_dist .* new_pos .+ (1 - relax_pos_dist) .* current_pos
            relax!(spin, parts, simulation, fraction_timestep, next_fraction_timestep)

            if ~has_intersection(collision)
                spin.position = new_pos
                found_solution = true
                break
            end
            relaxation = surface_relaxivity(simulation.geometry, simulation.properties, collision)
            if ~iszero(relaxation)
                transfer!.(spin.orientations, correct_for_timestep(relaxation, timestep))
            end

            permeability_prob = correct_for_timestep(permeability(simulation.geometry, simulation.properties, collision), timestep)
            passes_through = isone(permeability_prob) || !(iszero(permeability_prob) || rand() > permeability_prob)
            reflection = Reflection(collision, new_pos - current_pos, reflection.ratio_displaced, 
                reflection.time_moved + (1 - fraction_timestep) * use_distance * timestep, 
                reflection.distance_moved + norm(new_pos - current_pos) * use_distance, 
                passes_through
            )
            current_pos = spin.position = @. (current_pos * (1 - use_distance) + new_pos * use_distance)
            if ~isnothing(test_new_pos)
                push!(all_positions, current_pos)
            end

            sd = surface_density(simulation.geometry, simulation.properties, collision)
            if !iszero(sd) && rand() < stick_probability(sd, dwell_time(simulation.geometry, simulation.properties, collision), simulation.diffusivity, timestep)
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
        error("Bounced single particle for 10000 times in single step; terminating!")
    end
    if ~isnothing(test_new_pos)
        push!(all_positions, spin.position)
        return all_positions
    end
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
    spins::Vector{Spin{N}} = deepcopy.(snapshot.spins)

    times = propose_times(simulation, snapshot.time, new_time)

    # define next stopping times due to sequence, readout or times
    sequence_instants = Union{Nothing, InstantComponent}[next_instant(seq, current_time) for seq in simulation.sequences]
    sequence_times = MVector{N, Float64}([isnothing(i) ? Inf : get_time(i) for i in sequence_instants])

    for next_time in times
        # evolve all spins to next interesting time
        parts = SequencePart.(simulation.sequences, current_time, next_time)
        Threads.@threads for spin in spins
            draw_step!(spin, simulation, parts, next_time - current_time)
        end
        current_time = next_time

        # return final snapshot state
        if current_time != new_time && any(t -> t == current_time, sequence_times)
            components = SVector{N, Union{Nothing, InstantComponent}}([
                time == current_time ? instant : nothing 
                for (seq, instant, time) in zip(simulation.sequences, sequence_instants, sequence_times)
            ])
            apply!(components, spins)
            for (idx, ctime) in enumerate(sequence_times)
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    sequence_instants[idx] = next_instant(sequence, nextfloat(ctime))
                    sequence_times[idx] = get_time(sequence_instants[idx])
                end
            end
        end
    end
    return Snapshot(spins, current_time)
end

end