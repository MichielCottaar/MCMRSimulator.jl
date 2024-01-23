"""
Support for building NMR/MRI sequence diagrams.
"""
module Builder
include("building_blocks.jl")
include("sequence_builders.jl")
include("wait.jl")
include("gradients/gradients.jl")
include("pulses/pulses.jl")

import .BuildingBlocks: BuildingBlock, scanner_constraints!
export BuildingBlock, scanner_constraints!

import .Wait: WaitBlock
export WaitBlock

import .SequenceBuilders: SequenceBuilder, start_time, end_time, duration, TR
export SequenceBuilder, start_time, end_time, duration, TR

import .Gradients: PulsedGradient, InstantGradientBlock, qval, rise_time, flat_time, slew_rate, gradient_strength, bval
export PulsedGradient, InstantGradientBlock, qval, rise_time, flat_time, slew_rate, gradient_strength, bval

import .Pulses: InstantRFPulseBlock, flip_angle, phase
export InstantRFPulseBlock, flip_angle, phase

end