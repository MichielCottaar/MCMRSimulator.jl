@testset "Geometry plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"

@testset "3D geometry plots" begin

    for (name, geometry) in [
        ("bendy_cylinder", mr.BendyCylinder(control_point=[[0, 0, 0], [0, 0.3, 1]], radius=[0.3, 0.6], repeats=[2., 2., 2.], spline_order=2, closed=[0, 0, 1])),
        ("spheres", mr.Spheres(position=[[0, 0, 0], [0, 0, 1]], repeats=[3, 3, 3], radius=0.6)),
        ("cylinders", mr.Cylinders(position=[[0, 0], [1, 1]], repeats=[3, 3], radius=0.6)),
        ("walls", mr.Walls(position=0., repeats=3.)),
    ]
        @testset "Plot $name in 3D" begin
            function plot_mesh(fname)
                f = mr.plot_geometry(geometry)
                CairoMakie.save(fname, f.figure)
            end

            @visualtest plot_mesh "$dir/$(name)_3D.png" !isCI
        end
    end
end
@testset "Plot walls" begin
    geometry = [
        mr.Walls(repeats=1),
        mr.Walls(rotation=:y),
        mr.Walls(repeats=1, rotation=[1, 1, 1]),
    ]
    pp = mr.PlotPlane()

    function plot_walls(fname)
        f = mr.plot_geometry(pp, geometry)
        CairoMakie.save(fname, f.figure)
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
        f = mr.plot_off_resonance(pp, geometry)
        mr.plot_geometry!(f.axis, pp, geometry)
        CairoMakie.save(fname, f.figure)
    end

    @visualtest plot_myelinated_annuli "$dir/myelinated_annuli.png" !isCI
end

end
