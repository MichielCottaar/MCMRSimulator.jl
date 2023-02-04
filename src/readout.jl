"""
    Simulation(
        sequences; geometry=Obstruction[], diffusivity=0.,
        R1=0, T1=Inf, R2=0, T2=Inf, off_resonance=0, MT_fraction=0, permeability=0,, 
        timestep=nothing, gradient_precision=0.01, sample_displacement=5, sample_off_resonance=10,
    )

Defines the setup of the simulation and stores the output of the run.

# Argument
## General parameters:
- `sequences`: Vector of [`Sequence`](@ref) objects. During the spin random walk the simulation will keep track of the spin magnetisations for all of the provided sequences.
- `geometry`: Set of obstructions, which can be used to restrict the diffusion, produce off-resonance fields, alter the local T1/T2 relaxation, and as sources of magnetisation transfer.
- `diffusivity`: Rate of the random motion of the spins in um^2/ms.

## MRI properties
These parameters determine the evolution and relaxation of the spin magnetisation.
- `R1`/`T1`: sets the longitudinal relaxation rate (R1 in kHz) or relaxation time (T1=1/R1 in ms). This determines how fast the longitudinal magnetisation returns to its equilibrium value of 1.
- `R2`/`T2`: sets the transverse relaxation rate (R2 in kHz) or relaxation time (T2=1/R2 in ms). This determines how fast the transverse magnetisation is lost.
- `off_resonance`: Size of the off-resonance field in this voxel in kHz.
These MRI properties can be overriden for spins inside the [`BaseObstruction`](@ref) objects of the `geoemtry`.

## Collision parameters
These parameters determine how parameters behave when hitting the [`BaseObstruction`](@ref) objects of the `geoemtry`.
They can be overriden for individual [`BaseObstruction`] objects.
- `MT_fraction`: the fraction of magnetisation transfered between the obstruction and the water spin at each collision.
- `permeability`: the probability that the spin will pass through the obstruction.
Note that `MT_fraction` and `permeability` are internally adjusted to make their effect independent of the timestep (see [`correct_for_timestep`](@ref)).

## Timestep parameters
These parameters (`timestep`, `gradient_precision`, `sample_displacement`, `sample_off_resonance`) control the timepoints at which the simulation is evaluated.
The default values should work well.
For more details on how to adjust them, see [`TimeController`](@ref).

# Running the simulation
To run a [`Snapshot`](@ref) of spins through the simulations you can use one of the following functions:
- [`evolve`](@ref): evolves the spins in the snapshot until a single given time and returns that state in a new [`Snapshot`](@ref).
- [`trajectory`](@ref): returns full spin trajectory (recommended only for small number of spins).
- [`signal`](@ref): returns signal variation over time.
- [`readout`](@ref): returns the snapshots at the sequence readouts.
"""
struct Simulation{N, G<:Geometry, S<:Sequence}
    # N sequences, datatype T
    sequences :: SVector{N, S}
    diffusivity :: Float
    properties :: GlobalProperties
    geometry :: G
    time_controller::TimeController
    function Simulation(
        sequences, 
        diffusivity::Float,
        properties::GlobalProperties,
        geometry::Geometry,
        time_controller::TimeController
    )
        if isa(sequences, Sequence)
            sequences = [sequences]
        elseif length(sequences) == 0
            sequences = Sequence[]
        end
        nseq = length(sequences)
        geometry = Geometry(geometry)

        new{nseq, typeof(geometry), eltype(sequences)}(
            SVector{nseq}(sequences),
            diffusivity,
            properties,
            geometry,
            time_controller,
        )
    end
end

function Simulation(
    sequences;
    R1=NaN,
    T1=NaN,
    R2=NaN,
    T2=NaN,
    off_resonance=0.,
    MT_fraction=0.,
    permeability=0.,
    diffusivity=0.,
    geometry=Obstruction[],
    timestep=nothing,
    gradient_precision=0.01,
    sample_displacement=5,
    sample_off_resonance=10,
)
    geometry = Geometry(geometry)
    if isnothing(timestep)
        controller = TimeController(geometry; gradient_precision=gradient_precision, sample_displacement=sample_displacement, sample_off_resonance=sample_off_resonance)
    else
        controller = TimeController(timestep)
    end
    if iszero(diffusivity) && length(geometry) > 0
        @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
    end
    return Simulation(
        sequences, 
        diffusivity,
        GlobalProperties(; R1=R1, T1=T1, R2=R2, T2=T2, off_resonance=off_resonance, MT_fraction=MT_fraction, permeability=permeability),
        geometry,
        controller,
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

produces_off_resonance(sim::Simulation) = produces_off_resonance(sim.geometry)
propose_times(sim::Simulation, t_start, t_end) = propose_times(sim.time_controller, t_start, t_end, sim.sequences, sim.diffusivity)

"""
    readout(snapshot, simulation)

Evolves the spins in the [`Snapshot`](@ref) through the [`Simulation`](@ref).
Returns the [`Snapshot`](@ref) at every [`Readout`](@ref) in the simulated sequences during a single TR.
If no `TR` is explicitly selected, it will return the current TR if the snapshot has not passed any readouts and the next TR otherwise.
"""
function readout(spins, simulation::Simulation{N}) where {N}
    snapshot = _to_snapshot(spins, N)

    readout_times = Float[]
    per_seq_times = Vector{Float}[]

    for seq in simulation.sequences
        if length(seq.readout_times) == 0
            seq_times = Float[]
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
    for time in sort(readout_times)
        snapshot = evolve_to_time(snapshot, simulation, time)
        for index in 1:N
            if time in per_seq_times[index]
                push!(final_snapshots[index], get_sequence(snapshot, index))
            end
        end
    end
    return final_snapshots
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
        if get_time(snapshot) > times
            error("Trajectory final time is lower than the current time of the snapshot.")
        end
        times = propose_times(simulation, snapshot.time, times)
    end
    result = Array{typeof(snapshot)}(undef, size(times))
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))
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
        if get_time(snapshot) > times
            error("Trajectory final time is lower than the current time of the snapshot.")
        end
        times = propose_times(simulation, snapshot.time, times)
    end
    if N == 1
        result = Array{SpinOrientation}(undef, size(times))
    else
        result = Array{SVector{N, SpinOrientation}}(undef, size(times))
    end
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))
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
    evolve_to_time(snapshot, simulation, Float(new_time))
end