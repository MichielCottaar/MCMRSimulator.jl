module MRSimulator
import StaticArrays: SA_F64, MVector, SVector, @SVector
using LinearAlgebra
import Base
import RecipesBase: RecipesBase, @userplot, @recipe, @series
import CoordinateTransformations
import Rotations
import Colors


include("constants.jl")
include("spin.jl")
include("field.jl")
include("diffuse.jl")
include("sequence.jl")
include("readout.jl")
include("evolve.jl")
include("plot.jl")

export evolve, Microstructure, field, Sequence, RFPulse, Spin, vector, time, transverse, longitudinal, phase, Snapshot, evolve_iter
end

