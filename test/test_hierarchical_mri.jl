@testset "test_hierarchical_mri.jl" begin
    defaults = mr.GlobalProperties(R2=1.).mri
    for (obstruction, pos_in, pos_out) in [
        (mr.spheres(1., T1_inside=10., R2_inside=0.5), [0., 0, 0], [1, 2, 1.]),
        (mr.spheres(1., T1_inside=10., R2_inside=0.5, repeats=(5, 5, 5)), [5., 5, 0], [1, 2, 1.]),
        (mr.cylinders(1., T1_inside=10., R2_inside=0.5), [0., 0, 0], [1, 2, 1.]),
        (mr.annuli(0.5, 1., T1_inside=10., R2_inside=0.5), [0., 0, 0], [1, 2, 1.]),
    ]
        geom = mr.Geometry(obstruction)
        @testset "Test correct values in $(obstruction)" begin
            @test mr.isinside(geom, pos_in) > 0
            @test mr.isinside(geom, pos_out) == 0
            within = mr.inside_MRI_properties(geom, pos_in, defaults)
            @test mr.R1(within) == 0.1
            @test mr.T1(within) == 10.
            @test mr.R2(within) == 0.5
            @test mr.T2(within) == 2.
            @test iszero(mr.off_resonance(within))
            outside = mr.inside_MRI_properties(geom, pos_out, defaults)
            @test iszero(mr.R1(outside))
            @test isinf(mr.T1(outside))
            @test mr.R2(outside) == 1.
            @test mr.T2(outside) == 1.
            @test iszero(mr.off_resonance(outside))
        end
        @testset "Test simulation in $(obstruction)" begin
            sequence = mr.Sequence(pulses=[mr.InstantRFPulse(flip_angle=90., time=0.)], TR=100.)
            sim = mr.Simulation(sequence, R2=1., geometry=geom, diffusivity=1.)
            snap = mr.evolve([pos_in, pos_out], sim, 10.)
            within = snap[1].orientations[1]
            @test mr.transverse(within) ≈ exp(-5.)
            @test mr.longitudinal(within) ≈ 1 - exp(-1.)
            outside = snap[2].orientations[1]
            @test mr.transverse(outside) ≈ exp(-10.)
            @test abs(mr.longitudinal(outside)) < 1e-8
        end
    end
end
