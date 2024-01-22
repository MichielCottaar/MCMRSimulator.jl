module Gradients
include("pulsed_gradients.jl")

import .PulsedGradients: PulsedGradient, qval, rise_time, flat_time, slew_rate, gradient_strength
end