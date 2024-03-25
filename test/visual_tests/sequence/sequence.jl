@testset "Sequence plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"

@testset "Instantaneous gradients & pulses" begin
    function plot_perfect_dwi(fname)
        sequence = DWI(TE=80., bval=2., gradient=(type=:instant, ))
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_perfect_dwi "$dir/perfect_dwi.png" !isCI
end

@testset "Finite RF pulses" begin
    function plot_finite_rf(fname)
        times = -2:0.1:2
        raw_amp = sinc.(times)
        pulses = [
            GenericPulse(times .+ 2, raw_amp),
            InstantPulse(flip_angle=60.),
            GenericPulse(times .+ 2, raw_amp .* 2),
            InstantPulse(flip_angle=120.),
        ]

        sequence = build_sequence() do 
            Sequence(pulses) 
        end
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_finite_rf "$dir/finite_rf.png" !isCI
end

@testset "Finite gradients" begin
    function plot_finite_dwi(fname, single_gradient=false)
        sequence = DWI(bval=2., TE=80, TR=100, orientation=[0, -1, 1], scanner=Siemens_Prisma)
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence; single_gradient=single_gradient)
        CairoMakie.save(fname, f)
    end

    @visualtest fn->plot_finite_dwi(fn, true) "$dir/single_grad_dwi.png" !isCI
    @visualtest fn->plot_finite_dwi(fn, false) "$dir/multi_grad_dwi.png" !isCI
end

@testset "pulseseq example PRESS sequence" begin
    function plot_press(fname, single_gradient=false)
        fn = joinpath(dir, "..", "..", "example_pulseseq", "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
        sequence = mr.read_pulseq(fn)
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence; single_gradient=single_gradient)
        CairoMakie.save(fname, f)
    end

    @visualtest fn->plot_press(fn, true) "$dir/single_pulseq_press.png" !isCI
    @visualtest fn->plot_press(fn, false) "$dir/multi_pulseq_press.png" !isCI
end

end