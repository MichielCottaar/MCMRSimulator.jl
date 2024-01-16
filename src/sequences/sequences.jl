module Sequences
include("methods.jl")
include("instants.jl")
include("shapes.jl")
include("gradients.jl")
include("radio_frequency.jl")
include("main.jl")
include("pulseq.jl")
include("json.jl")

import .Methods: start_time, end_time, add_TR
import .Gradients: MRGradients, gradient, rotate_bvec
import .RadioFrequency: RFPulse, constant_pulse, amplitude, flip_angle, phase
import .Main: Sequence, SequencePart, effective_pulse
import .Main: previous_pulse, current_pulse, next_pulse
import .Main: previous_gradient, current_gradient, next_gradient
import .Main: previous_instant, current_instant, next_instant
import .PulseQ: read_pulseq
import .JSON: read_sequence, write_sequence
import .Instants: InstantComponent, Readout, InstantRFPulse, InstantGradient, apply!, qval, qvec
import .Shapes: control_points
end