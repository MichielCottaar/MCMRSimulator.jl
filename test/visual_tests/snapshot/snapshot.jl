@testset "Snapshot plots" begin
    isCI = get(ENV, "CI", "false") == "true"
    dir = @__DIR__
    @testset "Spins as dyads" begin
        function plot_dyads(fname; kind=:dyad, project=true)
            snapshot = mr.Snapshot(
                [
                    mr.Spin(position=[0, 0, 0], transverse=0.8, phase=0.)
                    mr.Spin(position=[1, 0, 1], transverse=0.8, phase=90.)
                    mr.Spin(position=[1, 1, 2], transverse=0.8, phase=180.)
                    mr.Spin(position=[0, 1, 3], transverse=0.8, phase=-90.)
                ]
            )
            f = Figure()
            if project
                plot_plane = mr.PlotPlane()
                plot(f[1, 1], plot_plane, snapshot; kind=kind)
            else
                plot(f[1, 1], snapshot; kind=kind)
            end
            CairoMakie.save(fname, f)
        end

        @visualtest plot_dyads "$dir/dyad_snapshot.png" !isCI
        @visualtest fn -> plot_dyads(fn, kind=:scatter) "$dir/scatter_snapshot.png" !isCI
        @visualtest fn -> plot_dyads(fn, project=false) "$dir/dyad_snapshot3D.png" !isCI
        @visualtest fn -> plot_dyads(fn, project=false, kind=:scatter) "$dir/scatter_snapshot3D.png" !isCI
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
            f = Figure()
            plot(f[1, 1], plot_plane, snapshot, kind=:image, ngrid=10)
            CairoMakie.save(fname, f)
        end

        @visualtest plot_image "$dir/image_snapshot.png" !isCI

    end
end
