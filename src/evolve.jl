function evolve_to_time(
    spin::Spin{N}, current_time::Float, new_time::Float,
    micro::Microstructure, timestep::Float=1., B0::Float=3.
) where {N}
    evolve_to_time(spin, current_time, new_time, micro, timestep, SVector{N, Float}(repeat([B0], N)))
end

"""
    evolve_to_time(spin, current_time, new_time, micro, timestep, B0)

Evolve a single spin to the next time of interest.
This takes into account both random diffusion of the spin's position
and relaxation of the MR spin orientation.
It is used internally when evolving [`Simulation`](@ref) objects.
"""
function evolve_to_time(
    spin::Spin{N}, current_time::Float, new_time::Float,
    micro::Microstructure, timestep::Float, B0::SVector{N, Float}
) where {N}
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return spin
    end
    index = div(current_time, timestep, RoundNearest)
    time_left = ((1//2 + index) * timestep) - current_time
    orient = MVector{N, SpinOrientation}(spin.orientations)

    function proc!(orient, pos, dt, micro=micro, B0=B0) 
        for idx in 1:N
            orient[idx] = relax(orient[idx], micro(pos), dt, B0[idx])
        end
    end

    if time_left > (new_time - current_time)
        # spin does not get to move at all
        proc!(orient, spin.position, new_time-current_time)
        return Spin(spin.position, SVector{N, SpinOrientation}(orient), spin.rng)
    end
    proc!(orient, spin.position, time_left)
    current_time = (index + 1//2) * timestep
    position::PosVector = spin.position

    # We need to move, so get random number generator from spin object
    old_rng_state = copy(Random.TaskLocalRNG())
    copy!(Random.TaskLocalRNG(), spin.rng)

    while new_time > (current_time + timestep)
        # Take full timesteps for a while
        if !isa(micro.diffusivity, ZeroField)
            position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
        end
        proc!(orient, position, timestep)
        current_time += timestep
    end

    if !isa(micro.diffusivity, ZeroField)
        position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
    end
    proc!(orient, position, new_time - current_time)

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

    B0_field = SVector{N, Float}(Float[s.B0 for s in simulation.sequences])
    nspins = length(spins)
    while true
        next = argmin(times)
        next_time::Float = times[next]

        # evolve all spins to next interesting time
        Threads.@threads for idx in 1:nspins
            spins[idx] = evolve_to_time(spins[idx], current_time, next_time, simulation.micro, simulation.timestep, B0_field)
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