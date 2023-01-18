function _relax_mult(spin::Spin{N}, dt, ct, sim::Simulation{N}) where {N}
    env = sim.micro(spin.position)

    function get_rng(sequence, orient)
        grad = get_gradient(spin.position, sequence, ct, dt + ct)
        relax(orient, env, dt, grad, sequence.scanner.B0)
    end

    orient = map(get_rng, sim.sequences, spin.orientations)
    Spin(spin.position, orient, spin.rng)
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

    if time_left > (new_time - current_time)
        # spin does not get to move at all
        return _relax_mult(spin, new_time-current_time, current_time, simulation)
    end
    spin = _relax_mult(spin, time_left, current_time, simulation)
    current_time = (index + 1//2) * simulation.timestep

    # We need to move, so get random number generator from spin object
    while new_time > (current_time + simulation.timestep)
        # Take full timesteps for a while
        if !isa(simulation.micro.diffusivity, ZeroField)
            spin = draw_step(spin, simulation.micro.diffusivity(spin.position), simulation.timestep, simulation.micro.geometry)
        end
        spin = _relax_mult(spin, simulation.timestep, current_time, simulation)
        current_time += simulation.timestep
    end

    if !isa(simulation.micro.diffusivity, ZeroField)
        spin = draw_step(spin, simulation.micro.diffusivity(spin.position), simulation.timestep, simulation.micro.geometry)
    end
   _relax_mult(spin, new_time - current_time, current_time, simulation)
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

    times = get_times(simulation, snapshot.time, new_time)
    # define next stopping times due to sequence, readout or times
    sequence_index = MVector{N, Int}(
        [next_pulse(seq, current_time) for seq in simulation.sequences]
    )
    sequence_times = MVector{N, Float}(
        [time(seq, index) for (seq, index) in zip(simulation.sequences, sequence_index)]
    )

    nspins = length(spins)
    for next_time in times
        # evolve all spins to next interesting time
        Threads.@threads for idx in 1:nspins
            spins[idx] = evolve_to_time(spins[idx], simulation, current_time, next_time)
        end
        current_time = next_time

        # return final snapshot state
        if current_time == new_time
            return Snapshot(spins, current_time)
        end

        # apply RF pulses
        if any(t -> t == current_time, sequence_times)
            components = SVector{N, Union{Nothing, SequenceComponent}}([
                time == current_time ? seq[index] : nothing 
                for (seq, index, time) in zip(simulation.sequences, sequence_index, sequence_times)
            ])
            spins = apply(components, spins)
            for (idx, ctime) in enumerate(times[2:end])
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    sequence_index[idx] += 1
                    sequence_times[idx] = time(sequence, sequence_index[idx])
                end
            end
        end
    end
end