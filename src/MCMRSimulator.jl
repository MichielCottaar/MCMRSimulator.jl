"""
This package supports the running of MR Monte Carlo simulations.

In these simulations hundreds of thousands or millions of [`Spin`](@ref) particles randomly diffuse through some tissue microstructure.
At each timepoint these spins are represented as a [`Snapshot`](@ref) object.
The spin diffusion is constrained by a set of [`Obstruction`](@ref) objects stored in a [`Geometry`](@ref) object.
The spins of these particles will be evolved based on the Bloch equations with the field strength and relaxation rates set by the local geometry
and the effect of one or more [`Sequence`](@ref) objects.
All these variables are combined into a single [`Simulation`](@ref) object. 
See [`Simulation`](@ref) for how to run the simulation.

Plotting support for the sequence and resulting signal is also available based on [Makie.jl](https://makie.juliaplots.org/stable/).
"""
module MCMRSimulator
import StaticArrays: SA_F64, MVector, SVector, @SVector, MMatrix, SA_F32, SMatrix
using LinearAlgebra
import Base
import Rotations
import CoordinateTransformations
import Colors
import Random
using UUIDs
import Distributions
import Accessors: @set
import Roots
import PlyIO
import Optim
import DataStructures: OrderedDict
using Statistics


include("constants.jl")
include("scanner.jl")
include("properties.jl")
include("geometry/geometry.jl")
include("spin.jl")
include("sequence/sequence.jl")
include("sequence_builder/sequence_builder.jl")
include("timestep.jl")
include("relax.jl")
include("readout.jl")
include("evolve.jl")
include("plot/plot.jl")

export Sequence, InstantRFPulse, Readout, InstantGradient, RFPulse
export MRGradients, gradient, rotate_bvec
export Scanner, B0, max_gradient, max_slew_rate
export Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export longitudinal, transverse, phase, vector, Spin, Snapshot, SpinOrientation, isinside, off_resonance, stuck, stuck_to
export TimeController, propose_times, size_scale
export Simulation, Microstructure, evolve, readout, signal, trajectory
export MRIProperties, CollisionProperties, GlobalProperties, ObstructionProperties, correct_for_timestep
export Annulus, annuli, Spiral, spirals, Cylinder, cylinders, Wall, walls, Sphere, spheres, Mesh, load_mesh
export random_positions_radii, TransformObstruction, Geometry, BoundingBox, get_rotation
export MT_fraction, permeability, surface_density, dwell_time, T1, R1, T2, R2, off_resonance
export previous_pulse, current_pulse, next_pulse, previous_gradient, current_gradient, next_gradient, previous_instant, current_instant, next_instant
export constant_pulse
export BuildingBlock, define_sequence, add_linear_diffusion_weighting
export read_pulseq, gradient_echo, spin_echo, dwi
export PlotPlane, plot_snapshot, image_snapshot, dyad_snapshot, plot_geometry, plot_trajectory2d, plot_trajectory3d, simulator_movie
end

