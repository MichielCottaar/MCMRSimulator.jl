module Gradients
include("integrate_gradients.jl")
include("pulsed_gradients.jl")

import .IntegrateGradients: qval, bval
import .PulsedGradients: PulsedGradient, rise_time, flat_time, slew_rate, gradient_strength
end