@testset "Sequence plots" begin

dir = @__DIR__

function plot_perfect_dwi(fname)
    sequence = mr.perfect_dwi(bval=2.)
    f = Figure()
    Axis(f[1, 1])
    plot!(sequence)
    CairoMakie.save(fname, f)
end

@visualtest plot_perfect_dwi "$dir/perfect_dwi.png"

end