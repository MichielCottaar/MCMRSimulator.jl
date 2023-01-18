"""
    Simulation(sequences, [microstructure]; R1=0., R2=0., off_resonance=0., diffusivity=0., geometry=Obstruction[], store_every=5., timestep=0.5)

Defines the setup of the simulation and stores the output of the run.

Note the use of units:
- all positions are in micrometers
- all times are in milliseconds
- magnetic fields are in Tesla

# Arguments
- `sequences::AbstractVector{Sequence}`: MR sequences to simulate. See [`Sequence`](@ref).
- `microstructure::Microstructure`: Description the tissue microstructure. A [`Microstructure`](@ref) object can be obtained directly or created using the `R1`, `R2`, `off_resonance`, `diffusivity`, and `geometry` flags.
- `timestep::Real`: Timestep of the spin random walk through the tissue in milliseconds (default: 0.5).

# Running the simulation
To run a [`Snapshot`](@ref) through the simulations you can use one of the following functions:
- [`evolve`](@ref): evolves the spins in the snapshot until a single given time and returns that state in a new [`Snapshot`](@ref).
- [`trajectory`](@ref): returns full spin trajectory (recommended only for small number of spins)
- [`signal`](@ref): returns signal variation over time
- [`readout`](@ref): returns the snapshots at the sequence readouts.

After creating the `Simulation` object it will not actually run.
To run the simulation for a given time, you can call:

    append!(simulation::Simulation, duration::Real)

This will run the simulation for `duration` milliseconds. 
During this time various intermediate states will be stored:
    - A [`Snapshot`](@ref)(`nsequences`) will be added to `simulation.regular` every `store_every` milliseconds.
    - A [`Snapshot`](@ref)(1) will be added to `simulation.readout[sequence_index]`, whenever a [`Readout`](@ref) object is encountered in the sequence.
    - The final state of the simulation will be added as a [`Snapshot`](@ref)(`nsequences`) object to `simulation.latest`, so that `simulation.latest[end]` always contains the final state of the total simulation.

`append!` can be called multiple times to continue the simulation further.
"""
struct Simulation{N, M<:Microstructure, S<:Sequence}
    # N sequences, datatype T
    sequences :: SVector{N, S}
    micro::M
    timestep::Float
    function Simulation(
        sequences, 
        micro::Microstructure; 
        timestep :: Real=0.5
    )
        if isa(sequences, Sequence)
            sequences = [sequences]
        elseif length(sequences) == 0
            sequences = Sequence[]
        end
        nseq = length(sequences)

        timestep = Float(timestep)

        new{nseq, typeof(micro), eltype(sequences)}(
            SVector{nseq}(sequences),
            micro,
            timestep,
        )
    end
end

function Simulation(
    sequences;
    R1=0.,
    R2=0.,
    diffusivity=0.,
    off_resonance=0.,
    geometry=Obstruction[],
    kwargs...
)
    return Simulation(
        sequences, Microstructure(R1=R1, R2=R2, diffusivity=diffusivity, off_resonance=off_resonance, geometry=geometry); kwargs...
    )
end

_to_snapshot(spins::Int, nseq) = _to_snapshot(Snapshot(spins), nseq)
_to_snapshot(spins::AbstractVector{<:Real}, nseq) = _to_snapshot(Spin(position=spins), nseq)
_to_snapshot(spins::AbstractVector{<:AbstractVector{<:Real}}, nseq) = _to_snapshot([Spin(position=pos) for pos in spins], nseq)
function _to_snapshot(spins::AbstractMatrix{<:Real}, nseq) 
    if size(spins, 2) != 3
        spins = transpose(spins)
    end
    @assert size(spins, 2) == 3
    _to_snapshot(Spin.(spins), nseq)
end
_to_snapshot(spins::Spin, nseq) = _to_snapshot([spins], nseq)
_to_snapshot(spins::AbstractVector{<:Spin}, nseq) = _to_snapshot(Snapshot(spins), nseq)
_to_snapshot(spins::Snapshot{1}, nseq::Int) = nseq == 1 ? spins : Snapshot(spins, nseq)
function _to_snapshot(spins::Snapshot{N}, nseq::Int) where {N}
    @assert nseq == N
    spins
end

produces_off_resonance(sim::Simulation) = produces_off_resonance(sim.micro)

"""
    readout(snapshot, simulation)

Evolves the spins in the [`Snapshot`](@ref) through the [`Simulation`](@ref).
Returns the [`Snapshot`](@ref) at every [`Readout`](@ref) in the simulated sequences during a single TR.
If no `TR` is explicitly selected, it will return the current TR if the snapshot has not passed any readouts and the next TR otherwise.
"""
function readout(spins, simulation::Simulation{N}) where {N}
    snapshot = _to_snapshot(spins, N)

    min_simulation_time = -1.

    get_times = Vector{Int}[]

    for seq in simulation.sequences
        current_TR = Int(div(time(snapshot), seq.TR, RoundDown))
        time_in_TR = time(snapshot) - current_TR * seq.TR
        times_readout = [time(r) for r in seq.pulses if isa(r, Readout)]
        if length(times_readout) == 0
            push!(get_times, Int[])
            continue
        end
        nskip = 0
        proposed_time = maximum(times_readout) + 1. + current_TR * seq.TR
        if minimum(times_readout) < time_in_TR
            proposed_time += seq.TR
            nskip = sum([t >= time_in_TR for t in times_readout])
        end
        get_times = push!(get_times, [(nskip + 1):(nskip + length(times_readout))...])
        if proposed_time > min_simulation_time
            min_simulation_time = proposed_time
        end
    end
    readouts = evolve_to_time(snapshot, simulation, Float(min_simulation_time))[2]
    return SVector{N}([
        r[indices] for (r, indices) in zip(readouts, get_times)
    ])
end

"""
    trajectory(snapshot, simulation, times=[TR])

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs at the requested times.
Returns a vector of [`Snapshot`](@ref) objects with the current state of each time in times.
When you are only interested in the signal at each timepoint, use [`signal`](@ref) instead.
"""
function trajectory(spins, simulation::Simulation{N}, times=nothing) where{N}
    snapshot = _to_snapshot(spins, N)
    if isnothing(times)
        times = simulation.sequences[1].TR
    end
    if isa(times, Real)
        if time(snapshot) > times
            error("Trajectory final time is lower than the current time of the snapshot.")
        end
        times = time(snapshot):simulation.timestep:times
    end
    result = Array{typeof(snapshot)}(undef, size(times))
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))[1]
        result[index] = snapshot
    end
    result
end

"""
    signal(snapshot, simulation, times=[TR])

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs the total signal at the requested times.
To get the full snapshot at each timepoint use [`trajectory`](@ref).
Returns a vector of [`SpinOrientation`](@ref) object with the total signal at each time in `times`.
"""
function signal(spins, simulation::Simulation{N}, times=nothing) where {N}
    snapshot = _to_snapshot(spins, N)
    if isnothing(times)
        times = simulation.sequences[1].TR
    end
    if isa(times, Real)
        if time(snapshot) > times
            error("Trajectory final time is lower than the current time of the snapshot.")
        end
        times = time(snapshot):simulation.timestep:times
    end
    if N == 1
        result = Array{SpinOrientation}(undef, size(times))
    else
        result = Array{SVector{N, SpinOrientation}}(undef, size(times))
    end
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))[1]
        if N == 1
            result[index] = SpinOrientation(snapshot)
        else
            result[index] = SVector{N}([SpinOrientation(get_sequence(snapshot, seq)) for seq in 1:N])
        end
    end
    result
end

"""
    evolve(snapshot, simulation[, new_time])

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) to a new time.
Returns a [`Snapshot`](@ref) at the new time, which can be used as a basis for further simulation.
By default it will simulate till the start of the next TR.
"""
function evolve(spins, simulation::Simulation{N}, new_time=nothing) where {N}
    snapshot = _to_snapshot(spins, N)
    if isnothing(new_time)
        TR = simulation.sequences[1].TR
        new_time = (div(nextfloat(snapshot.time), TR, RoundDown) + 1) * TR
    end
    evolve_to_time(snapshot, simulation, Float(new_time))[1]
end