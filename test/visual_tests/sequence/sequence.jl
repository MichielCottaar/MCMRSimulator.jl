@testset "Sequence plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"

@testset "Instantaneous gradients & pulses" begin
    function plot_perfect_dwi(fname)
        sequence = mr.dwi(bval=2., gradient_duration=0)
        f = mr.plot_sequence(sequence)
        CairoMakie.save(fname, f.figure)
    end

    @visualtest plot_perfect_dwi "$dir/perfect_dwi.png" !isCI
end

@testset "Finite RF pulses" begin
    function plot_finite_rf(fname)
        times = -2:0.1:2
        raw_amp = sinc.(times)
        pulses = [
            mr.RFPulse(times .+ 2, raw_amp),
            mr.InstantRFPulse(time=4., flip_angle=60.),
            mr.RFPulse(times .+ 6, raw_amp .* 2),
            mr.InstantRFPulse(time=8., flip_angle=120.),
        ]

        sequence = mr.Sequence(components=pulses, TR=30.)
        f = mr.plot_sequence(sequence)
        CairoMakie.save(fname, f.figure)
    end

    @visualtest plot_finite_rf "$dir/finite_rf.png" !isCI
end

@testset "Finite gradients" begin
    function plot_finite_dwi(fname)
        sequence = mr.dwi(bval=2., TE=80, TR=100, orientation=[0, -1, 1])
        f = mr.plot_sequence(sequence)
        CairoMakie.save(fname, f.figure)
    end

    @visualtest plot_finite_dwi "$dir/grad_dwi.png" !isCI
end

@testset "pulseseq example PRESS sequence" begin
    function plot_press(fname, single_gradient=false)
        fn = joinpath(dir, "..", "..", "example_pulseseq", "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
        sequence = mr.read_pulseq(fn)
        f = mr.plot_sequence(sequence; single_gradient=single_gradient)
        CairoMakie.save(fname, f.figure)
    end

    @visualtest fn->plot_press(fn, true) "$dir/single_pulseq_press.png" !isCI
    @visualtest fn->plot_press(fn, false) "$dir/multi_pulseq_press.png" !isCI
end

end