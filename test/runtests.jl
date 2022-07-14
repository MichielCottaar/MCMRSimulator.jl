using Test
import MRSimulator: Spin, ZeroField, Microstructure, evolve!, time, ConstantField
using StaticArrays

@testset "MRSimulator.jl" begin
    # Write your tests here.
    @testset "Simulate empty environment/sequence" begin
        spin = evolve!(Spin(), Microstructure(), 1.)
        @test spin.time == 1.
        @test spin.position == SA_F64[0., 0., 0.]
        @test spin.phase == 0.
    end
    @testset "Test constant off-resonance field" begin
        spin = evolve!(Spin(), Microstructure(off_resonance=ConstantField(2.)), 0.3)
        @test spin.time == 0.3
        @test spin.position == SA_F64[0., 0., 0.]
        @test spin.phase == 0.6
    end
end
