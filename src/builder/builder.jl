"""
Support for building NMR/MRI sequence diagrams.
"""
module Builder
include("building_blocks.jl")
include("sequence_builders.jl")
include("wait.jl")
include("gradients/gradients.jl")

import .BuildingBlocks: BuildingBlock, scanner_constraints!
export BuildingBlock, scanner_constraints!

import .Wait: WaitBlock
export WaitBlock

import .SequenceBuilders: SequenceBuilder, start_time, end_time, duration
export SequenceBuilder, start_time, end_time, duration

import .Gradients: PulsedGradient, qval, rise_time, flat_time, slew_rate, gradient_strength
export PulsedGradient, qval, rise_time, flat_time, slew_rate, gradient_strength

end