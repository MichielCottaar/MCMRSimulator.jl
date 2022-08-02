function evolve_to_time(spin::Spin, micro::Microstructure, current_time::Real, new_time::Real, B0=3.0)
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return spin
    end
    timestep = new_time - current_time
    if isa(micro.diffusivity, ZeroField)  # no diffusion
        return Spin(
            spin.position,
            relax(spin.orientation, micro(spin.position), timestep, B0)
        )
    elseif micro.geometry == Obstruction[]  # unrestricted diffusion
        step = draw_step(spin.position, micro.diffusivity(spin.position), timestep)
        return Spin(
            step.destination,
            relax(spin.orientation, micro(step), timestep, B0)
        )
    else
        steps = draw_step(spin.position, micro.diffusivity(spin.position), timestep, micro.geometry)
        orient = spin.orientation
        for step in steps
            orient = relax(orient, micro(step), timestep * step.timestep, B0)
        end
        return Spin(
            steps[end].destination,
            orient
        )
    end
end

function evolve_to_time(current::Snapshot, micro::Microstructure, new_time::Real, B0=3.0)
    if current.time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    nspins = length(current.spins)
    result = fill(Spin(), nspins)
    Threads.@threads for i = 1:nspins
        result[i] = evolve_to_time(current.spins[i], micro, current.time, new_time, B0)
    end
    Snapshot(result, new_time)
end

evolve_iter(snap::Snapshot, micro::Microstructure, sequence::Sequence; yield_every=1.0, B0=3.0) =
    Channel() do c
        sequence_index = 1
        readout_index = 1
        times = MVector{2,Float64}([
            time(sequence, sequence_index),
            (readout_index - 1) * yield_every
        ])
        while true
            next = argmin(times)
            snap = evolve_to_time(snap, micro, times[next], B0)
            if next == 1
                push!(c, snap)
                snap = Snapshot(
                    [apply(sequence[sequence_index], s) for s in snap.spins],
                    snap.time
                )
                push!(c, snap)
                sequence_index += 1
                times[1] = time(sequence, sequence_index)
            elseif next == 2
                push!(c, snap)
                readout_index += 1
                times[2] = (readout_index - 1) * yield_every
            end
        end
    end


function evolve_iter(spins::Vector{Spin}, micro::Microstructure, sequence::Sequence; yield_every=1.0, B0=3.0)
    evolve_iter(Snapshot(spins, 0.0), micro, sequence, yield_every=yield_every, B0=B0)
end

function evolve_iter(spin::Spin, micro::Microstructure, sequence::Sequence; yield_every=1.0, B0=3.0)
    evolve_iter([spin], micro, sequence, yield_every=yield_every, B0=B0)
end

function evolve(snap, micro::Microstructure, sequence::Sequence; yield_every=1.0, B0=3.0, nTR=1)
    max_time = nTR * sequence.TR
    all_snaps = Snapshot[]
    for snap in evolve_iter(snap, micro, sequence, yield_every=yield_every, B0=B0)
        if snap.time >= max_time && !(snap.time ≈ max_time)
            break
        end
        push!(all_snaps, snap)
        if snap.time ≈ max_time
            break
        end
    end
    all_snaps
end