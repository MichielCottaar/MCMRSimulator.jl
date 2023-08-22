@testset "Geometry plots" begin

dir = @__DIR__
isCI = get(ENV, "CI", "false") == "true"
@testset "Plot walls" begin
    geometry = [
        mr.walls(repeats=1),
        mr.walls(rotation=:y),
        mr.walls(repeats=1, rotation=[1, 1, 1]),
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
        mr.annuli(inner=0.2, outer=0.4, repeats=[1, 1], rotation=[0, 1, 1]),
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
        mr.annuli(inner=0.2, outer=0.4, repeats=[1, 1], rotation=[0, 1, 1], myelin=true)
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
