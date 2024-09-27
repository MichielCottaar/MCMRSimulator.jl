"""
This package supports the running of MR Monte Carlo simulations.

In these simulations hundreds of thousands or millions of [`Spin`](@ref) particles randomly diffuse through some tissue microstructure.
At each timepoint these spins are represented as a [`Snapshot`](@ref) object.
The spin diffusion is constrained by an [`ObstructionGroup`](@ref) (represented internally as [`FixedGeometry`](@ref))
The spins of these particles will be evolved based on the Bloch equations with the field strength and relaxation rates set by the local geometry
and the effect of one or more [`Sequence`](@ref) objects.
All these variables are combined into a single [`Simulation`](@ref) object. 
See [`Simulation`](@ref) for how to run the simulation.

Plotting support for the sequence and resulting signal is also available based on [Makie.jl](https://makie.juliaplots.org/stable/).
"""
module MCMRSimulator
include("constants.jl")
include("methods.jl")
include("properties.jl")
include("geometries/geometries.jl")
include("spins.jl")
include("timesteps.jl")
include("sequence_parts.jl")
include("simulations.jl")
include("relax.jl")
include("subsets.jl")
include("evolve.jl")
include("plot.jl")
include("cli/cli.jl")

import .Constants: gyromagnetic_ratio
export gyromagnetic_ratio

import .Methods: get_time, get_rotation
export get_time, get_rotation

import .Spins: position, longitudinal, transverse, phase, Spin, Snapshot, SpinOrientation, SpinOrientationSum, isinside, stuck, stuck_to, orientation
export position, longitudinal, transverse, phase, Spin, Snapshot, SpinOrientation, isinside, stuck, stuck_to, orientation

import .TimeSteps: TimeStep
export TimeStep

import .Simulations: Simulation, susceptibility_off_resonance
export Simulation, susceptibility_off_resonance

import .Subsets: Subset, get_subset
export Subset, get_subset

import .Evolve: evolve, readout
export evolve, readout

import .CLI: run_main, install_cli

import .Geometries: 
    ObstructionGroup, IndexedObstruction,
    Annulus, Annuli,
    Cylinder, Cylinders,
    Wall, Walls,
    Sphere, Spheres,
    Ring, BendyCylinder,
    Triangle, Mesh, nvolumes,
    load_mesh, fix, fix_susceptibility,
    random_positions_radii, write_geometry, read_geometry
export Annuli, Cylinders, Walls, Spheres, Mesh, load_mesh, fix, random_positions_radii, BendyCylinder

import .Geometries.Internal: BoundingBox, FixedGeometry, surface_relaxivity, surface_density, dwell_time, permeability
export BoundingBox

import .Properties: GlobalProperties, R1, R2, off_resonance, correct_for_timestep
export GlobalProperties, R1, R2, off_resonance

import .Plot: PlotPlane, plot_snapshot, plot_geometry, plot_trajectory, simulator_movie, plot_off_resonance
export PlotPlane, plot_snapshot, plot_geometry, plot_trajectory, simulator_movie, plot_off_resonance

import .Plot: plot_snapshot!, plot_geometry!, plot_trajectory!, plot_off_resonance!
export plot_snapshot!, plot_geometry!, plot_trajectory!, plot_off_resonance!


# Additional imports that will not be exported.
# They will not be of interest to most users.
import .Spins: FixedXoshiro
import .Methods: get_time, get_rotation, project
import .Evolve: draw_step!
import .Relax: relax!
end

