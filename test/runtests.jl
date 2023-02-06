using Test
import MCMRSimulator as mr
import MCMRSimulator: Float, SA
using StaticArrays
using LinearAlgebra
import Random
import SpecialFunctions: erf
using Statistics

all_tests = [
    "collisions",
    "evolve",
    "plots",
    "known_sequences",
    "meshes",
    "offresonance",
    "transfer",
    "permeability",
    "radio_frequency",
    "hierarchical_mri",
    "pulseseq",
    "various",
]

if length(ARGS) == 0
    tests = all_tests
elseif length(ARGS) == 1 && ARGS[1] == "no-plots"
    tests = symdiff(all_tests, ["plots"])
else
    tests = ARGS
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
