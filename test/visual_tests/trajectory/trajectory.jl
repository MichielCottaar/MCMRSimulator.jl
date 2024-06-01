@testset "Trajectory plots" begin
    isCI = get(ENV, "CI", "false") == "true"
    dir = @__DIR__
    @testset "3D plot" begin
        function plot_trajectory(fname)
            Random.seed!(1234)
            snapshot = mr.Snapshot([mr.Spin() for _ in 1:4])
            simulation = mr.Simulation(mr.dwi(bval=0.1, TE=80, gradient_duration=0), diffusivity=3.)
            trajectory = mr.readout(snapshot, simulation, 0:0.1:80, return_snapshot=true)
            f = Figure()
            mr.plot_trajectory3d(f[1, 1], trajectory, sequence=1, axis=(type=CairoMakie.Axis3,))
            CairoMakie.save(fname, f)
        end

        @visualtest plot_trajectory "$dir/trajectory_3d.png" !isCI
    end
    @testset "2D plot" begin
        function plot_trajectory(fname)
            Random.seed!(1234)
            snapshot = mr.Snapshot([mr.Spin() for _ in 1:4])
            simulation = mr.Simulation(mr.dwi(bval=0.1, TE=80, gradient_duration=0), diffusivity=3.)
            trajectory = mr.readout(snapshot, simulation, 0:0.1:80, return_snapshot=true)
            pp = mr.PlotPlane(sizex=30., sizey=30.)
            f = Figure()
            mr.plot_trajectory2d(f[1, 1], pp, trajectory, sequence=1)
            CairoMakie.save(fname, f)
        end

        @visualtest plot_trajectory "$dir/trajectory_2d.png" !isCI
    end
end