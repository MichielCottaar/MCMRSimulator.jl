"""
Support for building NMR/MRI sequence diagrams.
"""
module SequenceBuilder
include("building_blocks.jl")
include("gradients.jl")

import .BuildingBlocks: BuildingBlock, duration, to_components
import .Gradients: PulsedGradient
import .Gradients: rise_time, flat_time, duration, δ, qval, gradient_strength, slew_rate


end