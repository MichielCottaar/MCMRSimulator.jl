function evolve_to_time(
    spin::MultiSpin{N}, current_time::Real, new_time::Real,
    micro::Microstructure, timestep::Real=1., B0::Real=3.
) where {N}
    evolve_to_time(spin, current_time, new_time, micro, timestep, SVector{N}(repeat([B0], N)))
end

function evolve_to_time(
    spin::MultiSpin{N, T}, current_time::Real, new_time::Real,
    micro::Microstructure, timestep::Real, B0::SVector{N, <:Real}
) where {N, T}
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return spin
    end
    index = div(current_time, timestep, RoundNearest)
    time_left = ((0.5 + index) * timestep) - current_time
    orient = MVector{N, SpinOrientation{T}}(spin.orientations)

    function proc!(orient, pos, dt, micro=micro, B0=B0) 
        for idx in 1:N
            orient[idx] = relax(orient[idx], micro(pos), dt, B0[idx])
        end
    end

    if time_left > (new_time - current_time)
        # spin does not get to move at all
        proc!(orient, spin.position, new_time-current_time)
        return MultiSpin(spin.position, SVector{N, SpinOrientation{T}}(orient), spin.rng)
    end
    proc!(orient, spin.position, time_left)
    current_time = (index + 0.5) * timestep
    position = spin.position

    # We need to move, so get random number generator from spin object
    old_rng_state = copy(Random.TaskLocalRNG())
    copy!(Random.TaskLocalRNG(), spin.rng)

    while new_time > (current_time + timestep)
        # Take full timesteps for a while
        position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
        proc!(orient, position, timestep)
        current_time += timestep
    end

    position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
    proc!(orient, position, new_time - current_time)

    # Restore random number state
    final_rng_state = FixedXoshiro(copy(Random.TaskLocalRNG()))
    copy!(Random.TaskLocalRNG(), old_rng_state)
    MultiSpin(position, SVector{N, SpinOrientation{T}}(orient), final_rng_state)
end

function Base.append!(simulation::Simulation{N, T}, delta_time::T) where {N, T}
    spins = copy(simulation.latest[end].spins)
    current_time = simulation.latest[end].time
    new_time = current_time + delta_time
    sequence_index = MVector{N, Int}(
        [next_pulse(seq, current_time) for seq in simulation.sequences]
    )
    next_readout = div(current_time, simulation.store_every, RoundUp) * simulation.store_every
    times = [new_time, next_readout, time.(simulation.sequences, sequence_index)...]
    B0_field = SVector{N}([s.B0 for s in simulation.sequences])
    nspins = length(spins)
    while current_time < new_time
        next = argmin(times)
        next_time = times[next]
        for idx in 1:nspins
            spins[idx] = evolve_to_time(spins[idx], current_time, next_time, simulation.micro, simulation.timestep, B0_field)
        end
        current_time = next_time
        if times[1] == current_time
            push!(simulation.latest, MultiSnapshot(spins, current_time))
            break
        end
        if times[2] == current_time
            push!(simulation.regular, MultiSnapshot(copy(spins), current_time))
            times[2] += simulation.store_every
        end
        if any(t -> t == current_time, times[3:end])
            components = SVector{N, Union{Nothing, SequenceComponent}}([
                time == current_time ? seq[index] : nothing 
                for (seq, index, time) in zip(simulation.sequences, sequence_index, times[3:end])
            ])
            spins = apply(components, spins)
            for (idx, ctime) in enumerate(times[3:end])
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    if isa(sequence[sequence_index[idx]], Readout)
                        push!(simulation.readout[idx], Snapshot(get_sequence.(spins, idx), current_time))
                    end
                    sequence_index[idx] += 1
                    times[idx + 2] = time(sequence, sequence_index[idx])
                end
            end
        end
    end
    simulation
end