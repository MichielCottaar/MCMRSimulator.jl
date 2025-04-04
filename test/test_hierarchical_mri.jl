@testset "test_hierarchical_mri.jl" begin
    defaults = mr.GlobalProperties(R2=1., off_resonance=1.)
    for (obstruction, pos_in, pos_out) in [
        (mr.Spheres(radius=1., R1_inside=0.1, R2_inside=0.5, off_resonance_inside=0.1), [0., 0, 0], [1, 2, 1.]),
        (mr.Spheres(radius=1., R1_inside=0.1, R2_inside=0.5, repeats=(5, 5, 5), off_resonance_inside=0.1), [5., 5, 0], [1, 2, 1.]),
        (mr.Cylinders(radius=1., R1_inside=0.1, R2_inside=0.5, off_resonance_inside=0.1), [0., 0, 0], [1, 2, 1.]),
    ]
        @testset "Test correct values in $(typeof(obstruction))" begin
            @test mr.isinside(obstruction, pos_in) > 0
            @test mr.isinside(obstruction, pos_out) == 0
            @test mr.R1(pos_in, obstruction, defaults) == 0.1
            @test mr.R2(pos_in, obstruction, defaults) == 1.5
            @test mr.off_resonance(pos_in, obstruction, defaults) == (0., 1.1)
            @test iszero(mr.R1(pos_out, obstruction, defaults))
            @test mr.R2(pos_out, obstruction, defaults) == 1.
            @test mr.off_resonance(pos_out, obstruction, defaults) == (0., 1.)
        end
        @testset "Test simulation in $(typeof(obstruction))" begin
            sequence = GradientEcho(TE=100.)
            sim = mr.Simulation(sequence, R2=1., geometry=obstruction, diffusivity=1.)
            snap = mr.evolve([pos_in, pos_out], sim, 10.)
            within = snap[1].orientations[1]
            @test mr.transverse(within) ≈ exp(-15.)
            @test mr.longitudinal(within) ≈ 1 - exp(-1.)
            outside = snap[2].orientations[1]
            @test mr.transverse(outside) ≈ exp(-10.)
            @test abs(mr.longitudinal(outside)) < 1e-8
        end
    end
    @testset "Test correct values in annuli" begin
        geometry = mr.Annuli(inner=0.5, outer=1., R1_inner_volume=0.1, R1_outer_volume=0.2, R2_outer_volume=10.)
        inner = SVector{3}([0., 0., 0.])
        outer = SVector{3}([0.7, 0., 0.])
        outside = SVector{3}([1.2, 0., 0.])
        @test mr.R1(inner, geometry, defaults) == 0.1
        @test mr.R1(outer, geometry, defaults) == 0.2
        @test iszero(mr.R1(outside, geometry, defaults))

        @test mr.R2(inner, geometry, defaults) == 1.
        @test mr.R2(outer, geometry, defaults) == 11.
        @test mr.R2(outside, geometry, defaults) == 1.
    end
    @testset "Test length of inside geometry" begin
        geometry = mr.Annuli(inner=0.5, outer=1., R1_outer_volume=0.2)
        simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
        @test length(simulation.geometry) == 2
        @test length(simulation.inside_geometry) == 2

        geometry = mr.Annuli(inner=0.5, outer=1.)
        simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
        @test length(simulation.geometry) == 2
        @test length(simulation.inside_geometry) == 0

        geometry = mr.Annuli(inner=0.5, outer=1., R1_inner_volume=0.2)
        simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
        @test length(simulation.geometry) == 2
        @test length(simulation.inside_geometry) == 1
    end
end
