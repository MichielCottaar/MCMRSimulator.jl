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
include("scanners.jl")
include("properties.jl")
include("geometries/geometries.jl")
include("spin.jl")
include("sequences/sequences.jl")
include("sequence_builder/sequence_builder.jl")
include("timestep.jl")
include("relax.jl")
include("simulations.jl")
include("evolve.jl")
include("plot/plot.jl")

import .Constants: gyromagnetic_ratio
export gyromagnetic_ratio

import .Methods: get_time, get_rotation
export get_time, get_rotation

import .Sequences: Sequence, InstantRFPulse, Readout, InstantGradient, RFPulse, MRGradients, gradient, rotate_bvec, flip_angle
export Sequence, InstantRFPulse, Readout, InstantGradient, RFPulse, MRGradients, gradient, rotate_bvec, flip_angle

import .Scanners: Scanner, B0, max_gradient, max_slew_rate, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export Scanner, B0, max_gradient, max_slew_rate, Siemens_Connectom, Siemens_Prisma, Siemens_Terra

import .Spins: position, longitudinal, transverse, phase, Spin, Snapshot, SpinOrientation, isinside, stuck, stuck_to, orientation
export position, longitudinal, transverse, phase, Spin, Snapshot, SpinOrientation, isinside, stuck, stuck_to, orientation

import .Timestep: TimeController, propose_times
export TimeController, propose_times

import .Simulations: Simulation
export Simulation

import .Evolve: evolve, readout, signal, trajectory
export evolve, readout, signal, trajectory

import .Geometries: 
    ObstructionGroup, IndexedObstruction,
    Annulus, Annuli, annuli,
    Cylinder, Cylinders, cylinders,
    Wall, Walls, walls,
    Sphere, Spheres, spheres,
    Triangle, Mesh, mesh, split_mesh,
    load_mesh, fix, fix_susceptibility,
    random_positions_radii, write_geometry, read_geometry
export annuli, cylinders, walls, spheres, mesh, load_mesh, fix, random_positions_radii, split_mesh

import .Geometries.Internal: BoundingBox, FixedGeometry
export BoundingBox

import .Properties: GlobalProperties, R1, R2, off_resonance, permeability, surface_relaxivity, surface_density, dwell_time
export GlobalProperties, R1, R2, off_resonance, permeability, surface_relaxivity, surface_density, dwell_time

import .Sequences: previous_pulse, current_pulse, next_pulse, previous_gradient, current_gradient, next_gradient, previous_instant, current_instant, next_instant
export previous_pulse, current_pulse, next_pulse, previous_gradient, current_gradient, next_gradient, previous_instant, current_instant, next_instant

import .Sequences: constant_pulse, read_pulseq
export constant_pulse, read_pulseq

import .SequenceBuilder: BuildingBlock, define_sequence, add_linear_diffusion_weighting, gradient_echo, spin_echo, dwi
export BuildingBlock, define_sequence, add_linear_diffusion_weighting, gradient_echo, spin_echo, dwi

import .Plot: PlotPlane, plot_snapshot, image_snapshot, dyad_snapshot, plot_geometry, plot_trajectory2d, plot_trajectory3d, simulator_movie, plot_off_resonance
export PlotPlane, plot_snapshot, image_snapshot, dyad_snapshot, plot_geometry, plot_trajectory2d, plot_trajectory3d, simulator_movie, plot_off_resonance

import .Plot: plot_snapshot!, image_snapshot!, dyad_snapshot!, plot_geometry!, plot_trajectory2d!, plot_trajectory3d!, plot_off_resonance!
export plot_snapshot!, image_snapshot!, dyad_snapshot!, plot_geometry!, plot_trajectory2d!, plot_trajectory3d!, plot_off_resonance!


# Additional imports that will not be exported.
# They will not be of interest to most users.
import .Sequences: apply!, start_time, end_time, amplitude, SequencePart, qval, qvec
import .Spins: FixedXoshiro
import .Methods: get_time, get_rotation, phase, project
import .Evolve: draw_step!
import .Relax: relax!
end

