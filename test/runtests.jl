using Test
import MRSimulator: Spin, ZeroField, Microstructure, evolve!, time, ConstantField, GradientField
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
    @testset "Test gradient off-resonance field" begin
        micro = Microstructure(off_resonance=GradientField(SA_F64[1.5, 0., 0.], 2.))
        spin = evolve!(Spin(), micro, 0.3)
        @test spin.time == 0.3
        @test spin.position == SA_F64[0., 0., 0.]
        @test spin.phase == 0.6
        # Move spin and evolve a bit further in time
        spin.position = SA_F64[3., 0., 0.]
        spin = evolve!(spin, Microstructure(off_resonance=GradientField(SA_F64[1.5, 0., 0.], 2.)), 0.5)
        @test spin.time == 0.5
        @test spin.position == SA_F64[3., 0., 0.]
        @test spin.phase == 0.6 + (2. + 1.5 * 3.) * 0.2
    end
end
