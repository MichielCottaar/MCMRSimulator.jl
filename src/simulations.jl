module Simulations
import StaticArrays: SVector
import ..Geometries: ObstructionGroup, fix, fix_susceptibility
import ..Geometries.Internal: FixedGeometry, FixedObstruction, FixedSusceptibility
import ..Sequences: Sequence
import ..Timestep: TimeController, propose_times
import ..Spins: Spin, Snapshot, SpinOrientation
import ..Methods: get_time, B0
import ..Properties: GlobalProperties, R1, R2, off_resonance, correct_for_timestep

"""
    Simulation(
        sequences; geometry=[], diffusivity=0.,
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
These MRI properties can be overriden for spins inside the [`ObstructionGroup`](@ref) objects of the `geoemtry`.

## Collision parameters
These parameters determine how parameters behave when hitting the [`ObstructionGroup`](@ref) objects of the `geoemtry`.
They can be overriden for individual objects for each [`ObstructionGroup`].
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
- `evolve`: evolves the spins in the snapshot until a single given time and returns that state in a new [`Snapshot`](@ref).
- `readout`: evolves the spins to particular times in each TR and return the total signal at that time (or a [`Snapshot`](@ref)).
"""
struct Simulation{N, NG, G<:FixedGeometry{NG}, IG<:FixedGeometry, S<:Sequence, O<:FixedSusceptibility}
    # N sequences, datatype T
    sequences :: SVector{N, S}
    diffusivity :: Float64
    properties :: GlobalProperties
    geometry :: G
    inside_geometry :: IG
    susceptibility :: O
    time_controller::TimeController
    flatten::Bool
    function Simulation(
        sequences, 
        diffusivity::Float64,
        properties::GlobalProperties,
        geometry::FixedGeometry,
        inside_geometry::FixedGeometry,
        susceptibility::FixedSusceptibility,
        time_controller::TimeController,
        flatten::Bool,
    )
        nseq = length(sequences)

        new{nseq, length(geometry), typeof(geometry), typeof(inside_geometry), eltype(sequences), typeof(susceptibility)}(
            SVector{nseq}(sequences),
            diffusivity,
            properties,
            geometry,
            inside_geometry,
            susceptibility,
            time_controller,
            flatten
        )
    end
end

function Simulation(
    sequences;
    diffusivity=3,
    geometry=[],
    max_timestep=nothing,
    gradient_precision=1.,
    rf_rotation=1.,
    permeability=0.,
    surface_density=0.,
    dwell_time=0.,
    surface_relaxivity=0.,
    R1=0.,
    R2=0.,
    off_resonance=0.,
)
    flatten = false
    if isa(sequences, Sequence)
        sequences = [sequences]
        flatten = true
    elseif length(sequences) == 0
        sequences = Sequence[]
    end
    susceptibility = fix_susceptibility(geometry)
    geometry = fix(geometry; permeability=permeability, density=surface_density, dwell_time=dwell_time, relaxivity=surface_relaxivity)
    last_interesting_inside = findlast(g->~all(all(getproperty(g.volume, s) .== v) for (s, v) in ((:R1, 0), (:R2, 0), (:off_resonance, 0))), geometry)
    if isnothing(last_interesting_inside)
        inside_geometry = ()
    else
        inside_geometry = geometry[1:last_interesting_inside]
    end

    default_properties = GlobalProperties(; R1=R1, R2=R2, off_resonance=off_resonance)
    max_B0 = iszero(length(sequences)) ? 0 : maximum(B0.(sequences))
    controller = TimeController(geometry, susceptibility, max_B0, diffusivity; max_timestep=max_timestep, gradient_precision=gradient_precision, rf_rotation=rf_rotation)
    if iszero(diffusivity) && length(geometry) > 0
        @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
    end
    return Simulation(
        sequences, 
        Float64(diffusivity),
        default_properties,
        geometry,
        inside_geometry,
        susceptibility,
        controller,
        flatten,
    )
end

function Base.show(io::IO, sim::Simulation{N}) where {N}
    function print_geometry()
        print(io, "Geometry(")
        for obstruction in sim.geometry
            print(io, string(obstruction), ", ")
        end
        print(io, ")")
    end
    if get(io, :compact, false)
        seq_text = sim.flatten ? "single sequence" : "$N sequences"
        print(io, "Simulation($seq_text, ")
        print_geometry()
        print(io, ", D=$(sim.diffusivity)um^2/ms, $(sim.properties))")
    else
        print(io, "Simulation(")
        print_geometry()
        print(io, ", D=$(sim.diffusivity)um^2/ms, $(sim.properties)):\n")
        print(io, "$N sequences:\n")
        for seq in sim.sequences
            print(io, seq)
        end
    end
end

function Snapshot(nspins::Integer, simulation::Simulation{N}, bounding_box=500; kwargs...) where {N}
    Snapshot(nspins, bounding_box, simulation.geometry; nsequences=N, kwargs...)
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

for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(spin_or_pos, simulation::Simulation)
            $symbol(spin_or_pos, simulation.geometry, simulation.properties)
        end
    end
end

end