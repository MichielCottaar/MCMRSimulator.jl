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
    function test_func(seq, time, amplitude)
        grad = mr.gradient(seq, time)
        @test grad[1] ≈ amplitude[1] atol=1e-12
        @test grad[2] ≈ amplitude[2] atol=1e-12
        @test grad[3] ≈ amplitude[3] atol=1e-12
    end
    sequence = mr.Sequence(mr.BuildingBlock([grad_rot_def, grad_rot_x, grad_rot_y, grad_unrot_def, grad_unrot_x, grad_unrot_y]))

    test_func(sequence, 0.5, [0.5, 0., 0])
    test_func(sequence, 2., [1., 0., 0])
    test_func(sequence, 6., [1., 0., 0])
    test_func(sequence, 10., [0., 1., 0])
    test_func(sequence, 14., [1., 0., 0])
    test_func(sequence, 18., [1., 0., 0])
    test_func(sequence, 22., [0., 1., 0])

    rotated_y = mr.rotate_bvec(sequence, [0, 1., 0])
    test_func(rotated_y, 0.5, [0, 0.5, 0])
    test_func(rotated_y, 2., [0., 1., 0])
    test_func(rotated_y, 6., [0., 1., 0])
    test_func(rotated_y, 10., [-1., 0., 0])
    test_func(rotated_y, 14., [1., 0., 0])
    test_func(rotated_y, 18., [1., 0., 0])
    test_func(rotated_y, 22., [0., 1., 0])

    rotated_z = mr.rotate_bvec(sequence, [0, 0., 1])
    test_func(rotated_z, 0.5, [0., 0., 0.5])
    test_func(rotated_z, 2., [0., 0., 1])
    test_func(rotated_z, 6., [0., 0., 1])
    test_func(rotated_z, 10., [0., 1., 0])
    test_func(rotated_z, 14., [1., 0., 0])
    test_func(rotated_z, 18., [1., 0., 0])
    test_func(rotated_z, 22., [0., 1., 0])
end
end
