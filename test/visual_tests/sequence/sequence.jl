@testset "Sequence plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"

@testset "Instantaneous gradients & pulses" begin
    function plot_perfect_dwi(fname)
        sequence = mr.perfect_dwi(bval=2.)
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_perfect_dwi "$dir/perfect_dwi.png" !isCI
end

@testset "Finite gradients" begin
    function plot_finite_dwi(fname, single_gradient=false)
        sequence = mr.dwi(bval=2.)
        f = Figure()
        Axis(f[1, 1])
        plot!(sequence; single_gradient=single_gradient)
        CairoMakie.save(fname, f)
    end

    @visualtest fn->plot_perfect_dwi(fn, true) "$dir/single_grad_dwi.png" !isCI
    @visualtest fn->plot_perfect_dwi(fn, false) "$dir/multi_grad_dwi.png" !isCI
end

end