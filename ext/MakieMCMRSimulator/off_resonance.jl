module OffResonance
using Makie
import StaticArrays: SVector
import MCMRSimulator.Plot: PlotPlane, plot_off_resonance, plot_off_resonance!
import MCMRSimulator.Geometries.Internal: FixedSusceptibility, susceptibility_off_resonance
import MCMRSimulator.Geometries: fix_susceptibility

@Makie.recipe(Plot_Off_Resonance, plot_plane, geometry) do scene
    Makie.Theme(
        colormap=:viridis,
        ngrid=400
    )
end

function Makie.plot!(por::Plot_Off_Resonance)
    plot_plane = por[1]
    raw_geometry = por[2]
    susc = @lift raw_geometry isa FixedSusceptibility ? raw_geometry : fix_susceptibility($raw_geometry)

    dims = @lift -0.5:(1/$(por[:ngrid])):0.5
    xx_1d = @lift $dims * $plot_plane.sizex
    yy_1d = @lift $dims * $plot_plane.sizey
    pos_plane = @lift broadcast(
        (x, y) -> SVector{3}([x, y, 0.]),
        reshape($xx_1d, length($xx_1d), 1),
        reshape($yy_1d, 1, length($yy_1d)),
    )
    pos_orig = @lift inv($plot_plane.transformation).($pos_plane)
    field = @lift map(p->susceptibility_off_resonance($susc, p), $pos_orig)
    Makie.image!(por, xx_1d, yy_1d, field, colormap=por[:colormap])
    por
end

end