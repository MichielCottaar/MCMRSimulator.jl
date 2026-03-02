module Sequences

include("base.jl")
include("pulseq.jl")
include("sequence_io.jl")

import ..Base: Sequence, BuildingBlock, GradientWaveform, RFPulse, ADC, duration
export Sequence, BuildingBlock, GradientWaveform, RFPulse, ADC, duration

end