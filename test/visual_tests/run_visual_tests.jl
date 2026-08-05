@testset "Visual regression tests" begin
using CairoMakie
using VisualRegressionTests
using Gtk

include("geometry/geometry.jl")
include("snapshot/snapshot.jl")
include("trajectory/trajectory.jl")
end