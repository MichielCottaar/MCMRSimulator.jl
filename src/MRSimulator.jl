module MRSimulator
import StaticArrays: SA_F64, MVector, SVector, @SVector, MMatrix
using LinearAlgebra
import Base
import CoordinateTransformations
import Rotations
import Colors
import Random


include("constants.jl")
include("spin.jl")
include("field.jl")
include("diffuse.jl")
include("sequence/sequence.jl")
include("readout.jl")
include("evolve.jl")
include("plot.jl")

export Simulation, Microstructure, field, Sequence, RFPulse, Readout, InstantGradient, 
    longitudinal, transverse, phase, vector, Repeated, Transformed, Sphere, Cylinder, Wall,
    perfect_dwi
end

