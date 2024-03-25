using Test
import MCMRSimulator as mr
using StaticArrays
using LinearAlgebra
import Random
import SpecialFunctions: erf
using Statistics
using MRIBuilder

all_tests = [
    "collisions",
    "evolve",
    "known_sequences",
    "meshes",
    "offresonance",
    "transfer",
    "permeability",
    "radio_frequency",
    "hierarchical_mri",
    "sequence_IO",
    "various",
    "subsets",
    "plots",
    "cli",
]

if length(ARGS) == 0
    tests = all_tests
elseif length(ARGS) == 1 && ARGS[1] == "no-plots"
    tests = symdiff(all_tests, ["plots"])
else
    tests = ARGS
end


"""
    correct_collisions(movement, geometry)

Splits the given movement from point A to point B into multiple steps that bounce off the given obstructions.
This function assumes perfect reflection rather than the diffuse reflection used in [`draw_step`](@ref).
It is used to test the collision detection and resolution, but not actually used in the simulations.
"""
function correct_collisions(start, dest, geometry)
    simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
    spin = mr.Spin(nsequences=0, position=start)
    parts = SVector{0, mr.SequencePart}()
    return mr.draw_step!(spin, simulation, parts, 1., dest)
end

@testset "MCMRSimulator tests" begin
    for test in tests
        if test == "plots"
            include("visual_tests/run_visual_tests.jl")
        else
            include("test_$test.jl")
        end
    end
end
