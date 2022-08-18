"""
    Simulation(spins, sequences, [microstructure]; R1=0., R2=0., off_resonance=0., diffusivity=0., geometry=Obstruction[], store_every=5., timestep=0.5)

Defines the setup of the simulation and stores the output of the run.

Note the use of units:
- all positions are in micrometers
- all times are in milliseconds
- magnetic fields are in Tesla

# Arguments
- `spins`: One of [`Snapshot`](@ref). The latter allows all simulated sequences to start out in with a different spin orientation, while the former assumes them to be the same (note that the spin path is always the same across all sequences).
    - When starting from scratch the easiest way to create these is using one of the methods described in [`Snapshot`](@ref)
    - When starting from an existing simulation, you can access the final state of the simulation using `simulation.latest[end]`. 
        To reset the timer to zero, just extract the spins (`simulation.latest[end].spins`).
        To only extract the final result for a single sequence use [`get_sequence`](@ref) (e.g., `get_sequence(simulation.latest[end], sequence_index)`).
- `sequences::AbstractVector{Sequence}`: MR sequences to simulate. See [`Sequence`](@ref).
- `microstructure::Microstructure`: Description the tissue microstructure. A [`Microstructure`](@ref) object can be obtained directly or created using the `R1`, `R2`, `off_resonance`, `diffusivity`, and `geometry` flags.
- `store_every::Real`: Determines how often a snapshot should be stored (default: 5 milliseconds).
- `timestep::Real`: Timestep of the spin random walk through the tissue (default: 0.5 milliseconds).

# Running the simulation
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
struct Simulation{N, M<:Microstructure}
    # N sequences, datatype T
    sequences :: SVector{N, <:Sequence}
    micro::M
    timestep::Float
    regular :: AbstractVector{Snapshot{N}}
    store_every :: Float
    readout :: SVector{N, AbstractVector{Snapshot{1}}}
    latest :: Vector{Snapshot{N}}
    function Simulation(
        spins, 
        sequences, 
        micro::Microstructure; 
        store_every :: Real=5., 
        timestep :: Real=0.5
    )
        if isa(sequences, Sequence)
            sequences = [sequences]
        elseif length(sequences) == 0
            sequences = Sequence[]
        end
        nseq = length(sequences)

        store_every = Float(store_every)
        timestep = Float(timestep)

        if isa(spins, Spin)
            spins = [spins]
        end
        if isa(spins, AbstractVector{<:Spin})
            spins = Snapshot(spins)
        elseif isa(spins, AbstractVector{<:Spin})
            spins = Snapshot(spins)
        end
        if isa(spins, Snapshot)
            spins = Snapshot(spins, nseq)
        end

        new{nseq, typeof(micro)}(
            SVector{nseq}(sequences),
            micro,
            timestep,
            Snapshot{nseq}[],
            store_every,
            [Snapshot[] for _ in 1:nseq],
            [spins]
        )
    end
end

function Simulation(
    spins, 
    sequences;
    R1=0.,
    R2=0.,
    diffusivity=0.,
    off_resonance=0.,
    geometry=Obstruction[],
    kwargs...
)
    return Simulation(
        spins, sequences, 
        Microstructure(R1=R1, R2=R2, diffusivity=diffusivity, off_resonance=off_resonance, geometry=geometry); kwargs...
    )
end
