@testset "Visual regression tests" begin
using CairoMakie
using VisualRegressionTests
using Gtk

include("geometry/geometry.jl")
include("sequence/sequence.jl")
include("snapshot/snapshot.jl")
include("trajectory/trajectory.jl")
end