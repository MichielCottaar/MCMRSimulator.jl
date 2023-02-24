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

"""
    correct_collisions(movement, geometry)

Splits the given movement from point A to point B into multiple steps that bounce off the given obstructions.
This function assumes perfect reflection rather than the diffuse reflection used in [`draw_step`](@ref).
It is used to test the collision detection and resolution, but not actually used in the simulations.
"""
correct_collisions(to_try :: mr.Movement, geometry) = correct_collisions(to_try, mr.Geometry(geometry))
correct_collisions(to_try :: mr.Movement, geometry::mr.Geometry{0}) = [to_try]
function correct_collisions(to_try :: mr.Movement, geometry::mr.Geometry)
    steps = mr.Movement[]
    collision = mr.empty_collision
    while true
        collision = mr.detect_collision(to_try, geometry, collision)
        if collision === mr.empty_collision
            push!(steps, to_try)
            break
        end
        new_pos = collision.distance * to_try.destination + (1 - collision.distance) * to_try.origin
        push!(steps, mr.Movement(
            to_try.origin,
            new_pos,
            collision.distance * to_try.timestep
        ))
        if length(steps) >= 100
            error()
        end
        direction = to_try.destination .- to_try.origin
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction
        new_dest = new_pos .+ reflection / norm(reflection) * norm(direction) * (1 - collision.distance)
        to_try = mr.Movement(
            new_pos,
            new_dest,
            (1 - collision.distance) * to_try.timestep
        )
    end
    steps
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
