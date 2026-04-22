using Test
import MCMRSimulator as mr
using StaticArrays
using LinearAlgebra
import Random
import SpecialFunctions: erf
using Statistics
import Logging

Logging.disable_logging(Logging.Info)

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

It is used to test the collision detection and resolution (i.e., `Evolve.draw_step!`), but not actually used in the simulations.
"""
function correct_collisions(start, dest, geometry)
    simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
    spin = mr.Spin(nsequences=0, position=start)
    parts = mr.SequenceParts.MultSequencePart(1., mr.SequenceParts.SequencePart[], mr.SequenceParts.InstantSequencePart{SVector{0}}([]))
    B0s = SVector{0, Float64}()
    return mr.Evolve.draw_step!(spin, simulation, parts, B0s, SVector{3, Float64}(dest))
end

build_empty_sequence(duration) = mr.SequenceParts.SequenceWaveform((([], []), ([], []), ([], [])), [], [], [], 10.)

@testset "MCMRSimulator tests" begin
    for test in tests
        if test == "plots"
            include("visual_tests/run_visual_tests.jl")
        else
            include("test_$test.jl")
        end
    end
end
