"""
This package supports the running of MR Monte Carlo simulations.

In these simulations hundreds of thousands or millions of particles randomly diffuse through some tissue microstructure
constrained by a series of [`Obstruction`](@ref) objects.
The spins of these particles will be evolved by the effect of one or more [`Sequence`](@ref) objects as well as
spatially varying [`Field`](@ref) objecs describing the R1, R2, diffusivity, and off-resonance fields.
Off-resonance fields can also be generated by the [`Obstruction`](@ref) objects themselves.
All of these variables are combined into a single [`Simulation`](@ref) object. 
See [`Simulation`](@ref) for how to run the simulation.

Plotting support for the output is also available based on [Makie.jl](https://makie.juliaplots.org/stable/).
"""
module MRSimulator
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


include("constants.jl")
include("scanner.jl")
include("spin.jl")
include("field.jl")
include("geometry/geometry.jl")
include("sequence/sequence.jl")
include("microstructure.jl")
include("relax.jl")
include("readout.jl")
include("evolve.jl")
include("plot/plot.jl")

export Sequence, RFPulse, Readout, InstantGradient, create_gradients
export create_gradients, get_gradient, LinearGradients, StepWiseGradients, rotate_bvec
export dwi
export Scanner, Siemens_Connectom, Siemens_Prisma, Siemens_Terra
export longitudinal, transverse, phase, vector, Spin, Snapshot, SpinOrientation, isinside, off_resonance
export Simulation, Microstructure, evolve, readout, signal, trajectory, TransformObstruction
export Annulus, annuli, random_annuli, Spiral, spirals, random_spirals, Cylinder, cylinders, random_cylinders, Wall, walls, Sphere, spheres, random_spheres, Mesh, box_mesh
export PlotPlane, plot_snapshot, image_snapshot, dyad_snapshot, plot_geometry, plot_trajectory2d, plot_trajectory3d
end

