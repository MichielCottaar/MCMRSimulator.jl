@testset "test_gen_sequence.jl" begin

@testset "test rotating bvecs" begin
    times = [0, 1, 3, 4]
    amplitudes = [0, 1, 1, 0]
    grad_rot_def = mr.MRGradients(times, amplitudes, apply_bvec=true)
    grad_rot_x = mr.MRGradients(times, [[a, 0, 0] for a in amplitudes], apply_bvec=true)
    grad_rot_y = mr.MRGradients(times, [[0, a, 0] for a in amplitudes], apply_bvec=true)
    grad_unrot_def = mr.MRGradients(times, amplitudes)
    grad_unrot_x = mr.MRGradients(times, [[a, 0, 0] for a in amplitudes])
    grad_unrot_y = mr.MRGradients(times, [[0, a, 0] for a in amplitudes])
    sequence = mr.Sequence(mr.BuildingBlock([grad_rot_def, grad_rot_x, grad_rot_y, grad_unrot_def, grad_unrot_x, grad_unrot_y]))
    @test mr.amplitude(sequence.gradient, 0.5) == [0.5, 0., 0]
    @test mr.amplitude(sequence.gradient, 2.) == [1., 0., 0]
    @test mr.amplitude(sequence.gradient, 6.) == [1., 0., 0]
    @test mr.amplitude(sequence.gradient, 10.) == [0., 1., 0]
    @test mr.amplitude(sequence.gradient, 14.) == [1., 0., 0]
    @test mr.amplitude(sequence.gradient, 18.) == [1., 0., 0]
    @test mr.amplitude(sequence.gradient, 22.) == [0., 1., 0]

    rotated_y = mr.apply_bvec(sequence, [0, 1., 0])
    @test mr.amplitude(rotated_y.gradient, 0.5) == [0, 0.5, 0]
    @test mr.amplitude(rotated_y.gradient, 2.) == [0., 1., 0]
    @test mr.amplitude(rotated_y.gradient, 6.) == [0., 1., 0]
    @test mr.amplitude(rotated_y.gradient, 10.) == [-1., 0., 0]
    @test mr.amplitude(rotated_y.gradient, 14.) == [1., 0., 0]
    @test mr.amplitude(rotated_y.gradient, 18.) == [1., 0., 0]
    @test mr.amplitude(rotated_y.gradient, 22.) == [0., 1., 0]

    rotated_z = mr.apply_bvec(sequence, [0, 0., 1])
    @test mr.amplitude(rotated_z.gradient, 0.5) == [0., 0., 0.5]
    @test mr.amplitude(rotated_z.gradient, 2.) == [0., 0., 1]
    @test mr.amplitude(rotated_z.gradient, 6.) == [0., 0., 1]
    @test mr.amplitude(rotated_z.gradient, 10.) == [0., 1., 0]
    @test mr.amplitude(rotated_z.gradient, 14.) == [1., 0., 0]
    @test mr.amplitude(rotated_z.gradient, 18.) == [1., 0., 0]
    @test mr.amplitude(rotated_z.gradient, 22.) == [0., 1., 0]
end
end
