@testset "Sequence plots" begin
    isCI = get(ENV, "CI", "false") == "true"
    dir = @__DIR__
    pulseq_dir = joinpath(@__DIR__, "..", "..", "pulseq")

    for (name, filename) in [
        ("instant_diffusion", "dwi_te_80_bval_2_instant_diffusion_time_40.seq"),
        ("delta_0", "dwi_te_80_bval_0.3_gradient_δ_0_Δ_nothing.seq"),
        ("regular_dwi", "dwi_te_80_bval_2.seq"),
        ("PRESS", "07_PRESS_offcenter.seq"),
    ]
        @testset "$name" begin
            function plot_single_sequence(fname)
                sequence = mr.read_pulseq(joinpath(pulseq_dir, filename))
                f = Figure()
                plot(f[1, 1], sequence)
                save_rgb(fname, f)
            end

            @visualtest plot_single_sequence "$dir/$(name).png" !isCI
        end
    end
end
