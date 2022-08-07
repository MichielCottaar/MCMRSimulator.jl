module MRSimulator
import StaticArrays: SA_F64, MVector, SVector, @SVector
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

export evolve, Microstructure, field, Sequence, RFPulse, Spin, vector, time, transverse, longitudinal, phase, Snapshot, evolve_iter
end

