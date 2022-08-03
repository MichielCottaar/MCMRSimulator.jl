function evolve_to_time(
    orient::SpinOrientation, trajectory::StepTrajectory,
    current_time::Real, new_time::Real,
    micro::Microstructure, B0::Real=3.
)
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return orient
    end
    index = Int(div(current_time, trajectory.timestep, RoundNearest)) + 1
    time = (index - 1) * trajectory.timestep
    next = minimum((time + 0.5 * trajectory.timestep, new_time))
    dt = next - current_time
    orient = relax(orient, micro(trajectory[index]), dt, B0)
    current_time = next

    while new_time > current_time
        index += 1
        next = minimum((current_time + trajectory.timestep, new_time))
        dt = next - current_time
        orient = relax(orient, micro(trajectory[index]), dt, B0)
        current_time = next
    end
    orient
end


function evolve_TR(spin::Spin, readout::SpinReadout, micro::Microstructure, trajectory::StepTrajectory)
    sequence = readout.parent.sequence
    sequence_index = 1
    readout_index = 1
    times = MVector{2,Float64}([
        time(sequence, sequence_index),
        time(readout.parent.all[readout_index])
    ])
    orient = spin.orientation
    current_time = 0.
    while readout_index <= length(readout.parent.all)
        next = argmin(times)
        orient = evolve_to_time(orient, trajectory, current_time, times[next], micro, sequence.B0)
        current_time = times[next]
        if next == 1
            orient = apply(sequence[sequence_index], orient)
            sequence_index += 1
            times[1] = time(sequence, sequence_index)
        elseif next == 2
            readout[readout_index] = Spin(trajectory[current_time], orient)
            readout_index += 1
            if readout_index <= length(readout.parent.all)
                times[2] = readout.parent.all[readout_index].time
            end
        end
    end
end

function evolve_TR(
    spin::Spin, 
    readouts::AbstractVector{SpinReadout},
    micro::Microstructure,
    timestep=0.1,
)
    trajectory = StepTrajectory(spin.position, timestep, micro.diffusivity, micro.geometry)
    [evolve_TR(spin, readout, micro, trajectory) for readout in readouts]
end


function evolve_TR(snap :: Snapshot, sequences::AbstractVector{Sequence}, micro::Microstructure; store_every=1.0)
    readouts = [
        TR_Readout(s, length(snap), store_every)
        for s in sequences
    ]
    for (index, spin) in enumerate(snap.spins)
        evolve_TR(spin, [SpinReadout(r, index) for r in readouts], micro)
    end
    readouts
end

function evolve_TR(snap :: Snapshot, sequence::Sequence, micro::Microstructure; kwargs...)
    evolve_TR(snap, [sequence], micro; kwargs...)[1]
end

function evolve_TR(spins :: AbstractVector{Spin}, sequences, micro::Microstructure; kwargs...)
    evolve_TR(Snapshot(spins, 0.), sequences, micro; kwargs...)
end

function evolve_TR(spin :: Spin, sequences, micro::Microstructure; kwargs...)
    evolve_TR([spin], sequences, micro; kwargs...)
end