module Sequences

include("base.jl")
include("pulseq_io/pulseq_io.jl")
include("pulseq.jl")
include("sequence_io.jl")

import .Base: Sequence, BaseBuildingBlock, BuildingBlock, GradientWaveform, RFPulse, ADC, duration, InstantGradient, InstantPulse
import .SequenceIO: read_sequence, write_sequence

end