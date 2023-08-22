module Sequences
include("gradient_echo.jl")
include("spin_echo.jl")
import .GradientEcho: gradient_echo
import .SpinEcho: spin_echo, dwi
end