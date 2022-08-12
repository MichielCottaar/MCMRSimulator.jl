"""
This package supports the running of MR Monte Carlo simulations.

In these simulations hundreds of thousands or millions of particles randomly diffuse through some tissue microstructure
constrained by a series of [`Obstruction`](@ref) objects.
The spins of these particles will be evolved by the effect of one or more [`Sequence`](@ref) objects as well as
spatially varying [`Field`](@ref) objecs describing the R1, R2, diffusivity, and off-resonance fields.
All of these variables are combined into a single [`Simulation`](@ref) object, which can then be evolved using [`append!`](@ref)(simulation).
Plotting support for the output is also available based on [Makie.jl](https://makie.juliaplots.org/stable/) (see [`plot!`](@ref)).
"""
module MRSimulator
import StaticArrays: SA_F64, MVector, SVector, @SVector, MMatrix
using LinearAlgebra
import Base
import CoordinateTransformations
import Rotations
import Colors
import Random
using UUIDs


include("constants.jl")
include("spin.jl")
include("field.jl")
include("geometry/geometry.jl")
include("sequence/sequence.jl")
include("microstructure.jl")
include("readout.jl")
include("evolve.jl")
include("plot.jl")

export Simulation, Microstructure, field, Sequence, RFPulse, Readout, InstantGradient, 
    longitudinal, transverse, phase, vector, Repeated, Transformed, Sphere, Cylinder, Wall,
    perfect_dwi
end

