@testset "Visual regression tests" begin
using CairoMakie
using VisualRegressionTests
using Test
using Gtk
isCI = get(ENV, "CI", "false") == "true"

include("sequence/sequence.jl")
end