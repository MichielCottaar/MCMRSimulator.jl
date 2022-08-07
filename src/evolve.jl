function evolve_to_time(
    spin::MultiSpin{N}, current_time::Real, new_time::Real,
    micro::Microstructure, timestep::Real=1., B0::Real=3.
) where {N}
    evolve_to_time(spin, current_time, new_time, micro, timestep, SVector{N}(repeat([B0], N)))
end

function evolve_to_time(
    spin::MultiSpin{N}, current_time::Real, new_time::Real,
    micro::Microstructure, timestep::Real, B0::SVector{N, <:Real}
) where {N}
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    proc(orient, pos, dt) = SVector{N}([relax(o, micro(pos), dt, b) for (o, b) in zip(orient, B0)])
    if new_time == current_time
        return spin
    end
    index = div(current_time, timestep, RoundNearest)
    time_left = ((0.5 + index) * timestep) - current_time
    if time_left > (new_time - current_time)
        # spin does not get to move at all
        orient = proc(spin.orientations, spin.position, new_time-current_time)
        return MultiSpin(
            spin.position,
            orient,
            spin.rng
        )
    end
    orient = proc(spin.orientations, spin.position, time_left)
    current_time = (index + 0.5) * timestep
    position = spin.position

    # We need to move, so get random number generator from spin object
    old_rng_state = copy(Random.TaskLocalRNG())
    copy!(Random.TaskLocalRNG(), spin.rng)

    while new_time > (current_time + timestep)
        # Take full timesteps for a while
        position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
        orient = proc(orient, position, timestep)
        current_time += timestep
    end

    position = draw_step(position, micro.diffusivity(position), timestep, micro.geometry)
    orient = proc(orient, position, new_time - current_time)

    # Restore random number state
    final_rng_state = FixedXoshiro(copy(Random.TaskLocalRNG()))
    copy!(Random.TaskLocalRNG(), old_rng_state)
    MultiSpin(position, orient, final_rng_state)
end

function Base.append!(simulation::Simulation{N, M}, delta_time::Real) where {N, M}
    snap = simulation.latest
    new_time = snap.time + delta_time
    sequence_index = MVector{N, Int}(
        [next_pulse(seq, snap.time) for seq in simulation.sequences]
    )
    next_readout = div(snap.time, simulation.store_every, RoundUp) * simulation.store_every
    times = [new_time, next_readout, time.(simulation.sequences, sequence_index)...]
    B0_field = SVector{N}([s.B0 for s in simulation.sequences])
    spins = snap.spins
    while snap.time < new_time
        next = argmin(times)
        snap = MultiSnapshot(
            SVector{M}([evolve_to_time(s, snap.time, times[next], simulation.micro, simulation.timestep, B0_field) for s in snap.spins]),
            times[next]
        )
        simulation.latest = snap
        if times[1] == snap.time
            break
        end
        if times[2] == snap.time
            push!(simulation.regular, snap)
            times[2] += simulation.store_every
        end
        if any(t -> t == snap.time, times[3:end])
            components = SVector{N, Union{Nothing, SequenceComponent}}([
                time == snap.time ? seq[index] : nothing 
                for (seq, index, time) in zip(simulation.sequences, sequence_index, times[3:end])
            ])
            snap = apply(components, snap)
            for (idx, ctime) in enumerate(times[3:end])
                if ctime == snap.time
                    sequence = simulation.sequences[idx]
                    if isa(sequence[sequence_index[idx]], Readout)
                        push!(simulation.readout[idx], get_sequence(snap, idx))
                    end
                    sequence_index[idx] += 1
                    times[idx + 2] = time(sequence, sequence_index[idx])
                end
            end
        end
    end
    simulation
end