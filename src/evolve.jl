"""
Defines the functions that run the actual simulation:
- [`readout`](@ref): get total signal or [`Snapshot`](@ref) at any [`Readout`](@ref) objects in the sequences.
- [`evolve`](@ref): Return a single [`Snapshot`](@ref) with the state of the simulation at a given time. This snapshot can be used as initialisation for further runs.

All of these functions call [`evolve_to_time`](@ref) under the hood to actually run the simulation.
"""
module Evolve
import StaticArrays: SVector, MVector
import LinearAlgebra: norm, ⋅
import MRIBuilder: Sequence, readout_times, TR
import Rotations
import ..SequenceParts: SequencePart
import ..Methods: get_time
import ..Spins: @spin_rng, Spin, Snapshot, stuck, SpinOrientationSum, get_sequence, orientation, SpinOrientation
import ..Simulations: Simulation, _to_snapshot
import ..Relax: relax!, transfer!
import ..Properties: GlobalProperties, correct_for_timestep, stick_probability
import ..Subsets: Subset, get_subset
import ..Geometries.Internal: Reflection, detect_intersection, empty_intersection, has_intersection, surface_relaxivity, permeability, surface_density, direction, previous_hit, dwell_time, empty_reflection

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
function readout(spins, simulation::Simulation{N}, new_readout_times=nothing; bounding_box=500, skip_TR=0, nTR=nothing, return_snapshot=false, noflatten=false, subset=Subset()) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)

    if isnothing(new_readout_times)
        if iszero(N)
            nreadout_per_TR = 0
        else
            nreadout_per_TR = maximum((length(readout_times(seq)) for seq in simulation.sequences))
        end
    else
        nreadout_per_TR = length(new_readout_times)
    end
    actual_readout_times = Set{Float64}()
    use_nTR = isnothing(nTR) ? 1 : nTR
    single_subset = ~(subset isa AbstractVector)
    subsets_vector = single_subset ? [subset] : subset
    store_times = fill(-1., (N, nreadout_per_TR, use_nTR))

    for (i, seq) in enumerate(simulation.sequences)
        current_TR = Int(div(get_time(snapshot), TR(seq), RoundDown))
        time_in_TR = get_time(snapshot) - current_TR * TR(seq)

        rt = isnothing(new_readout_times) ? readout_times(seq) : new_readout_times
        if rt isa Number
            rt = [rt]
        end
        if iszero(length(rt))
            continue
        end

        if iszero(skip_TR) && minimum(rt) < time_in_TR
            skip_TR = 1
        end

        for (j, relative_time) in enumerate(rt)
            for k in 1:use_nTR
                new_readout_time = ((skip_TR + k - 1) * TR(seq)) + relative_time
                if new_readout_time < snapshot.time
                    warn("Skipping readouts scheduled before the current snapshot. Did you provide negative `new_readout_times`?")
                    continue
                end
                push!(actual_readout_times, new_readout_time)
                store_times[i, j, k] = new_readout_time
            end
        end
    end

    return_type = return_snapshot ? Snapshot{1} : SpinOrientationSum
    result = convert(Array{Union{Nothing, return_type}}, fill(nothing, (size(store_times)..., length(subsets_vector))))

    for time in sort!([actual_readout_times...])
        snapshot = evolve_to_time(snapshot, simulation, time)
        for index in eachindex(IndexCartesian(), store_times)
            if time != store_times[index]
                continue
            end
            single_snapshot = get_sequence(snapshot, index[1])
            for (index_selected, select) in enumerate(subsets_vector)
                selected = get_subset(single_snapshot, simulation, select)
                value = return_snapshot ? selected : SpinOrientationSum(selected)
                result[index, index_selected] = value
            end
        end
    end

    # Only include Nothing in the return type if there are any nothings in the data.
    if ~any(isnothing.(result))
        result = convert(Array{return_type}, result)
    end

    if noflatten
        return result
    end

    return result[
        # remove first dimension if simulation was initalised with scalar sequence
        simulation.flatten ? 1 : :,
        # remove second dimension if there are no sequences with multiple Readouts or new_readout_times is set to a scalar number
        (isnothing(new_readout_times) && isone(nreadout_per_TR)) || readout_times isa Number ? 1 : :,
        # remove third dimension if nTR is not explicitly set by the user
        isnothing(nTR) ? 1 : :,
        single_subset ? 1 : :,
    ]
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
        first_TR = TR(simulation.sequences[1])
        if !all(TR(s) ≈ first_TR for s in simulation.sequences)
            error("Cannot evolve snapshot for a single TR, because the simulation contains sequences with different TRs. Please set a `new_time` explicitly.")
        end
        new_time = (div(nextfloat(snapshot.time), first_TR, RoundDown) + 1) * first_TR
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
function draw_step!(spin :: Spin{N}, simulation::Simulation{N}, parts::SVector{N, SequencePart}, timestep :: Float64, B0s::SVector{N, Float64}, test_new_pos=nothing) where {N}
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
    fraction_timestep = 0.

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
                td = dwell_time(simulation.geometry, spin.reflection)
                fraction_stuck = -log(rand()) * td / timestep
                relax!(spin, parts, simulation, fraction_timestep, min(one(Float64), fraction_timestep + fraction_stuck), B0s)
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
            next_fraction_timestep = fraction_timestep + (1 - fraction_timestep) * use_distance

            # spin relaxation
            relax_pos_dist = rand() * use_distance
            spin.position = map((p1, p2) -> relax_pos_dist * p1 + (1 - relax_pos_dist) * p2, new_pos, current_pos)
            relax!(spin, parts, simulation, fraction_timestep, next_fraction_timestep, B0s)

            if ~has_intersection(collision)
                spin.position = new_pos
                found_solution = true
                break
            end
            relaxation = surface_relaxivity(simulation.geometry, collision)
            if ~iszero(relaxation)
                transfer!.(spin.orientations, correct_for_timestep(relaxation, timestep))
            end

            permeability_prob = correct_for_timestep(permeability(simulation.geometry, collision), timestep)
            passes_through = isone(permeability_prob) || !(iszero(permeability_prob) || rand() > permeability_prob)
            reflection = Reflection(collision, new_pos - current_pos, reflection.ratio_displaced, 
                reflection.time_moved + (1 - fraction_timestep) * use_distance * timestep, 
                reflection.distance_moved + norm(new_pos - current_pos) * use_distance, 
                passes_through
            )
            current_pos = spin.position = map((p1, p2) -> use_distance * p1 + (1 - use_distance) * p2, new_pos, current_pos)
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
    if new_time == current_time
        return snapshot
    end
    spins::Vector{Spin{N}} = deepcopy.(snapshot.spins)

    linear_sequences = LinearSequence(simulation.sequences, current_time, new_time; max_timestep=simulation.max_timestep)
    B0s = map(seq->seq.scanner.B0, simulation.sequences)
    if iszero(N)
        if isinf(simulation.max_timestep)
            Nparts = 1
        else
            Nparts = Int(div(new_time - current_time, simulation.max_timestep, RoundUp))
        end
    else
        Nparts = length(linear_sequences[1].finite_parts)
    end

    for index_part in 1:Nparts
        # apply instant components
        for index_seq in 1:N
            if !isnothing(linear_sequences[index_seq].instant_pulses[index_part])
                # apply pulse to all spins
                pulse = linear_sequences[index_seq].instant_pulses[index_part]
                rotation = Rotations.RotationVec(
                    deg2rad(pulse.flip_angle) * cosd(pulse.phase),
                    deg2rad(pulse.flip_angle) * sind(pulse.phase),
                    0.
                )
                for spin in spins
                    orient = spin.orientations[index_seq]
                    new_orient = SpinOrientation(rotation * orientation(orient))
                    orient.longitudinal = new_orient.longitudinal
                    orient.transverse = new_orient.transverse
                    orient.phase = new_orient.phase
                end
            end
            if !isnothing(linear_sequences[index_seq].instant_gradient[index_part])
                # apply gradient to all spins
                grad = linear_sequences[index_seq].instant_gradient[index_part].qvec
                for spin in spins
                    new_phase = rad2deg(spin.position ⋅ grad)
                    spin.orientations[index_seq].phase += new_phase
                end
            end
        end

        parts = map(seq -> seq.finite_parts[index_part], linear_sequences)
        timestep = iszero(N) ? (new_time - current_time) / Nparts : parts[1].duration
        # evolve all spins to next interesting time
        Threads.@threads for spin in spins
            draw_step!(spin, simulation, parts, timestep, B0s)
        end
    end
    return Snapshot(spins, new_time)
end

end