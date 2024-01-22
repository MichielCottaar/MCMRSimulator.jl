"""
Support for building NMR/MRI sequence diagrams.
"""
module Builder
include("building_blocks.jl")
include("wait.jl")
include("sequence_builders.jl")
include("gradients/gradients.jl")

import .BuildingBlocks: BuildingBlock, scanner_constraints!
import .Wait: WaitBlock
import .SequenceBuilders: SequenceBuilder
import .Gradients: PulsedGradient

export BuildingBlock, scanner_constraints!
export PulsedGradient, WaitBlock, SequenceBuilder

end