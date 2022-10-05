function _proc!(orient::MVector{N}, pos, dt, ct, sim::Simulation{N}) where {N}
    for idx in 1:N
        sequence = sim.sequences[idx]
        grad = get_gradient(pos, sequence, ct, dt + ct)
        orient[idx] = relax(orient[idx], sim.micro(pos), dt, grad, sequence.scanner.B0)
    end
end

"""
    evolve_to_time(spin, simulation, current_time, new_time)

Evolve a single spin to the next time of interest.
This takes into account both random diffusion of the spin's position
and relaxation of the MR spin orientation.
It is used internally when evolving [`Simulation`](@ref) objects.
"""
function evolve_to_time(
    spin::Spin{N}, simulation::Simulation{N}, current_time::Float, new_time::Float
) where {N}
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return spin
    end
    index = Int(div(current_time, simulation.timestep, RoundNearest))
    time_left = ((1//2 + index) * simulation.timestep) - current_time
    orient = MVector{N, SpinOrientation}(spin.orientations)

    if time_left > (new_time - current_time)
        # spin does not get to move at all
        _proc!(orient, spin.position, new_time-current_time, current_time, simulation)
        return Spin(spin.position, SVector{N, SpinOrientation}(orient), spin.rng)
    end
    _proc!(orient, spin.position, time_left, current_time, simulation)
    current_time = (index + 1//2) * simulation.timestep
    position::PosVector = spin.position

    # We need to move, so get random number generator from spin object
    old_rng_state = copy(Random.TaskLocalRNG())
    copy!(Random.TaskLocalRNG(), spin.rng)

    while new_time > (current_time + simulation.timestep)
        # Take full timesteps for a while
        if !isa(simulation.micro.diffusivity, ZeroField)
            position = draw_step(position, simulation.micro.diffusivity(position), simulation.timestep, simulation.micro.geometry)
        end
        _proc!(orient, position, simulation.timestep, current_time, simulation)
        current_time += simulation.timestep
    end

    if !isa(simulation.micro.diffusivity, ZeroField)
        position = draw_step(position, simulation.micro.diffusivity(position), simulation.timestep, simulation.micro.geometry)
    end
    _proc!(orient, position, new_time - current_time, current_time, simulation)

    # Restore random number state
    final_rng_state = FixedXoshiro(copy(Random.TaskLocalRNG()))
    copy!(Random.TaskLocalRNG(), old_rng_state)
    Spin(position, SVector{N, SpinOrientation}(orient), final_rng_state)
end

"""
    evolve_to_time(snapshot, simulation, new_time)

Evolves the full [`Snapshot`](@ref) through the [`Simulation`](@ref) to the given `new_time`.
Multi-threading is used to evolve multiple spins in parallel.
This is used internally when calling any of the snapshot evolution methods (e.g., [`evolve`](@ref)).
"""
function evolve_to_time(snapshot::Snapshot{N}, simulation::Simulation{N}, new_time::Float) where {N}
    current_time::Float = snapshot.time
    if new_time < current_time
        error("New requested time ($(new_time)) is less than current time ($(snapshot.time)). Simulator does not work backwards in time.")
    end
    spins::Vector{Spin{N}} = copy(snapshot.spins)

    # define next stopping times due to sequence, readout or times
    sequence_index = MVector{N, Int}(
        [next_pulse(seq, current_time) for seq in simulation.sequences]
    )
    times = [new_time, time.(simulation.sequences, sequence_index)...]

    readouts = SVector{N, Vector{Snapshot{1}}}([Snapshot{1}[] for _ in 1:N])

    nspins = length(spins)
    while true
        next = argmin(times)
        next_time::Float = times[next]

        # evolve all spins to next interesting time
        Threads.@threads for idx in 1:nspins
            spins[idx] = evolve_to_time(spins[idx], simulation, current_time, next_time)
        end
        current_time = next_time

        # take care of any readouts
        if any(t -> t == current_time, times[2:end])
            for (idx, ctime) in enumerate(times[2:end])
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    if isa(sequence[sequence_index[idx]], Readout)
                        push!(readouts[idx], Snapshot(get_sequence.(spins, idx), current_time))
                        sequence_index[idx] += 1
                        times[idx + 1] = time(sequence, sequence_index[idx])
                    end
                end
            end
        end

        # return final snapshot state
        if times[1] == current_time
            return (Snapshot(spins, current_time), readouts)
        end

        # apply RF pulses
        if any(t -> t == current_time, times[2:end])
            components = SVector{N, Union{Nothing, SequenceComponent}}([
                time == current_time ? seq[index] : nothing 
                for (seq, index, time) in zip(simulation.sequences, sequence_index, times[2:end])
            ])
            spins = apply(components, spins)
            for (idx, ctime) in enumerate(times[2:end])
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    sequence_index[idx] += 1
                    times[idx + 1] = time(sequence, sequence_index[idx])
                end
            end
        end
    end
end