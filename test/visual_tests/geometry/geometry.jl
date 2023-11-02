@testset "Geometry plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"

@testset "Plot 3-dimensional mesh" begin
    geometry = mr.BendyCylinder(control_point=[[0, 0, 0], [0, 0.3, 1]], radius=[0.3, 0.6], repeats=[2., 2., 2.], spline_order=2, closed=[0, 0, 1])

    function plot_mesh(fname)
        f = Figure()
        mr.plot_geometry3d(f[1, 1], geometry)
        CairoMakie.save(fname, f)
    end

    function plot_mesh2(fname)
        f = Figure()
        mr.plot(f[1, 1], geometry)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_mesh "$dir/bendy_cylinder.png" !isCI
    @visualtest plot_mesh "$dir/bendy_cylinder.png" !isCI
end
@testset "Plot walls" begin
    geometry = [
        mr.Walls(repeats=1),
        mr.Walls(rotation=:y),
        mr.Walls(repeats=1, rotation=[1, 1, 1]),
    ]
    pp = mr.PlotPlane()

    function plot_walls(fname)
        f = Figure()
        mr.plot_geometry(f[1, 1], pp, geometry)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_walls "$dir/walls.png" !isCI
end
@testset "Plot annuli" begin
    geometry = [
        mr.Annuli(inner=0.2, outer=0.4, repeats=[1, 1], rotation=[0, 1, 1]),
    ]
    pp = mr.PlotPlane()

    function plot_annuli(fname)
        f = Figure()
        mr.plot_geometry(f[1, 1], pp, geometry)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_annuli "$dir/annuli.png" !isCI
end
@testset "Plot myelinated annuli" begin
    geometry = [
        mr.Annuli(inner=0.2, outer=0.4, repeats=[1, 1], rotation=[0, 1, 1], myelin=true)
    ]
    pp = mr.PlotPlane()

    function plot_myelinated_annuli(fname)
        f = Figure()
        mr.plot_off_resonance(f[1, 1], pp, geometry)
        mr.plot_geometry!(f[1, 1], pp, geometry)
        CairoMakie.save(fname, f)
    end

    @visualtest plot_myelinated_annuli "$dir/myelinated_annuli.png" !isCI
end

end
