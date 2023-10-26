module SequenceBuilder
include("building_blocks.jl")
include("define_sequence.jl")
include("diffusion.jl")
include("sequences/sequences.jl")

import .BuildingBlocks: BuildingBlock
import .DefineSequence: define_sequence
import .Diffusion: add_linear_diffusion_weighting, trapezium_gradient
import .Sequences: gradient_echo, spin_echo, dwi
end