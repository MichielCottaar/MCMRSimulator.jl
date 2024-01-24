module Pulses
include("properties.jl")
include("instant_pulses.jl")
include("constant_pulses.jl")

import .Properties: flip_angle, phase, amplitude, frequency
import .InstantPulses: InstantRFPulseBlock
import .ConstantPulses: ConstantPulse

end