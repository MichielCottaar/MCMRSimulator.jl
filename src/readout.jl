"""
    Simulation(
        sequences; geometry=Obstruction[], diffusivity=0.,
        R1=0, T1=Inf, R2=0, T2=Inf, off_resonance=0, MT_fraction=0, permeability=0,, 
        max_timestep=<geometry-based default>, gradient_precision=1, rf_rotation=1.,
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
- `surface_density`: Density of spins stuck on the surface relative to the volume density of hte free water.
- `dwell_time`: Typical time that spins will stay at the surface after getting stuck.
Note that `MT_fraction` and `permeability` are internally adjusted to make their effect independent of the timestep (see [`correct_for_timestep`](@ref)).

## Timestep parameters
These parameters (`max_timestep`, `gradient_precision`, and `rf_rotation`) control the timepoints at which the simulation is evaluated.
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
    flatten::Bool
    function Simulation(
        sequences, 
        diffusivity::Float,
        properties::GlobalProperties,
        geometry::Geometry,
        time_controller::TimeController,
        flatten::Bool,
    )
        nseq = length(sequences)

        new{nseq, typeof(geometry), eltype(sequences)}(
            SVector{nseq}(sequences),
            diffusivity,
            properties,
            geometry,
            time_controller,
            flatten
        )
    end
end

function Simulation(
    sequences;
    R1=NaN,
    T1=NaN,
    R2=NaN,
    T2=NaN,
    off_resonance=0,
    MT_fraction=0,
    permeability=0,
    diffusivity=0,
    surface_density=0,
    dwell_time=NaN,
    geometry=Obstruction[],
    max_timestep=nothing,
    gradient_precision=1.,
    rf_rotation=1.,
)
    flatten = false
    if isa(sequences, Sequence)
        sequences = [sequences]
        flatten = true
    elseif length(sequences) == 0
        sequences = Sequence[]
    end
    geometry = Geometry(geometry)
    default_properties = GlobalProperties(; 
        R1=R1, T1=T1, R2=R2, T2=T2, off_resonance=off_resonance, 
        MT_fraction=MT_fraction, permeability=permeability, surface_density=surface_density, dwell_time=dwell_time
    )
    max_B0 = iszero(length(sequences)) ? 0 : maximum(B0.(sequences))
    controller = TimeController(geometry, max_B0, diffusivity, default_properties; max_timestep=max_timestep, gradient_precision=gradient_precision, rf_rotation=rf_rotation)
    if iszero(diffusivity) && length(geometry) > 0
        @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
    end
    return Simulation(
        sequences, 
        Float(diffusivity),
        default_properties,
        geometry,
        controller,
        flatten,
    )
end

function Base.show(io::IO, sim::Simulation{N}) where {N}
    if get(io, :compact, false)
        seq_text = sim.flatten ? "single sequence" : "$N sequences"
        print(io, "Simulation($seq_text, $(sim.geometry), D=$(sim.diffusivity)um^2/ms, $(sim.properties))")
    else
        print(io, "Simulation($(sim.geometry), D=$(sim.diffusivity)um^2/ms, $(sim.properties)):\n")
        print(io, "$N sequences:\n")
        for seq in sim.sequences
            print(io, seq)
        end
    end
end

function Snapshot(nspins::Integer, simulation::Simulation, bounding_box=500; kwargs...)
    Snapshot(nspins, bounding_box, simulation.geometry, surface_density(simulation.properties); kwargs...)
end
_to_snapshot(spins::Int, simulation::Simulation, bounding_box) = _to_snapshot(Snapshot(spins, simulation, bounding_box), simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:Real}, simulation::Simulation, bounding_box) = _to_snapshot(Spin(position=spins), simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:AbstractVector{<:Real}}, simulation::Simulation, bounding_box) = _to_snapshot([Spin(position=pos) for pos in spins], simulation, bounding_box)
function _to_snapshot(spins::AbstractMatrix{<:Real}, simulation::Simulation, bounding_box) 
    if size(spins, 2) != 3
        spins = transpose(spins)
    end
    @assert size(spins, 2) == 3
    _to_snapshot([spins[i, :] for i in 1:size(spins, 1)], simulation, bounding_box)
end
_to_snapshot(spins::Spin, simulation::Simulation, bounding_box) = _to_snapshot([spins], simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:Spin}, simulation::Simulation, bounding_box) = _to_snapshot(Snapshot(spins), simulation, bounding_box)
_to_snapshot(spins::Snapshot{1}, simulation::Simulation{nseq}, bounding_box) where {nseq} = nseq == 1 ? spins : Snapshot(spins, nseq)
_to_snapshot(spins::Snapshot{N}, simulation::Simulation{N}, bounding_box) where {N} = spins

produces_off_resonance(sim::Simulation) = produces_off_resonance(sim.geometry)
propose_times(sim::Simulation, t_start, t_end) = propose_times(sim.time_controller, t_start, t_end, sim.sequences, sim.diffusivity)

"""
    readout(snapshot, simulation; bounding_box=<1x1x1 mm box>)

Evolves the spins in the [`Snapshot`](@ref) through the [`Simulation`](@ref).
Returns the [`Snapshot`](@ref) at every [`Readout`](@ref) in the simulated sequences during a single TR.
If no `TR` is explicitly selected, it will return the current TR if the snapshot has not passed any readouts and the next TR otherwise.

The return object depends on whether the simulation was created with a single sequence object or with a vector of sequences.
- For a single sequence object a vector of [`Snapshot`](@ref) objects is returned with a single snapshot for each [`Readout`](@ref) in the sequence.
- For a vector of sequences a vector of vectors of [`Snapshot`](@ref) objects is returned. Each element in the outer vector contains the result for a single sequence.
"""
function readout(spins, simulation::Simulation{N}; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)

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
    for time in sort(unique(readout_times))
        snapshot = evolve_to_time(snapshot, simulation, time)
        for index in 1:N
            if time in per_seq_times[index]
                push!(final_snapshots[index], get_sequence(snapshot, index))
            end
        end
    end
    if simulation.flatten
        return final_snapshots[1]
    else
        return final_snapshots
    end
end

"""
    trajectory(snapshot, simulation, times=[TR]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs at the requested times.
Returns a vector of [`Snapshot`](@ref) objects with the current state of each time in times.
When you are only interested in the signal at each timepoint, use [`signal`](@ref) instead.
"""
function trajectory(spins, simulation::Simulation{N}, times; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    result = Array{typeof(snapshot)}(undef, size(times))
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))
        result[index] = snapshot
    end
    result
end

"""
    signal(snapshot, simulation, times=[TR]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) and outputs the total signal at the requested times.
To get the full snapshot at each timepoint use [`trajectory`](@ref).

The return object depends on whether the simulation was created with a single sequence object or with a vector of sequences.
- For a single sequence object a vector of [`SpinOrientation`](@ref) object with the total signal at each time is returned.
- For a vector of sequences the return type is a (Nt, Ns) matrix of [`SpinOrientation`](@ref) objects for `Nt` timepoints and `Ns` sequences.
"""
function signal(spins, simulation::Simulation{N}, times; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    if simulation.flatten
        result = Array{SpinOrientation}(undef, size(times))
    else
        result = Array{SpinOrientation}(undef, (length(times), N))
    end
    for index in sortperm(times)
        snapshot = evolve_to_time(snapshot, simulation, Float(times[index]))
        if simulation.flatten
            result[index] = SpinOrientation(snapshot)
        else
            for seq in 1:N
                result[index, seq] = SpinOrientation(get_sequence(snapshot, seq))
            end
        end
    end
    result
end

"""
    evolve(snapshot, simulation[, new_time]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) to a new time.
Returns a [`Snapshot`](@ref) at the new time, which can be used as a basis for further simulation.
By default it will simulate till the start of the next TR.
"""
function evolve(spins, simulation::Simulation{N}, new_time=nothing; bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    if isnothing(new_time)
        TR = simulation.sequences[1].TR
        new_time = (div(nextfloat(snapshot.time), TR, RoundDown) + 1) * TR
    end
    evolve_to_time(snapshot, simulation, Float(new_time))
end