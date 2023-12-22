@testset "Snapshot plots" begin
    isCI = get(ENV, "CI", "false") == "true"
    dir = @__DIR__
    @testset "Spins as dyads or scatter plot" begin
        function plot_dyads(fname, kind)
            snapshot = mr.Snapshot(
                [
                    mr.Spin(position=[0, 0, 0], transverse=0.8, phase=0.)
                    mr.Spin(position=[1, 0, 1], transverse=0.8, phase=90.)
                    mr.Spin(position=[1, 1, 2], transverse=0.8, phase=180.)
                    mr.Spin(position=[0, 1, 3], transverse=0.8, phase=-90.)
                ]
            )
            plot_plane = mr.PlotPlane()
            f = mr.plot_snapshot(plot_plane, snapshot; dyadlength=1., kind=kind)
            CairoMakie.save(fname, f.figure)
        end

        @visualtest fn -> plot_dyads(fn, :dyad) "$dir/dyad_snapshot.png" !isCI
        @visualtest fn -> plot_dyads(fn, :scatter) "$dir/scatter_snapshot.png" !isCI
    end
    @testset "Spins as image" begin
        function plot_image(fname)
            Random.seed!(1234)
            positions = rand(SVector{3, Float64}, 30000)
            snapshot = mr.Snapshot(
                [
                    mr.Spin(position=pos, transverse=pos[1], phase=pos[2] * 360.)
                    for pos in positions
                ]
            )
            plot_plane = mr.PlotPlane()
            f = mr.plot_snapshot(plot_plane, snapshot, ngrid=10, kind=:image)
            CairoMakie.save(fname, f.figure)
        end

        @visualtest plot_image "$dir/image_snapshot.png" !isCI

    end
end
