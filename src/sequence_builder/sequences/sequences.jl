module Sequences
include("gradient_echo.jl")
include("spin_echo.jl")
include("stimulated_echo.jl")
include("SSFP.jl")

import .GradientEcho: gradient_echo
import .SpinEcho: spin_echo, dwi
import .SSFP: dwssfp
import .StimulatedEcho: stimulated_echo, dwste
end
