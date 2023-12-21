module OffResonance
using Makie
import StaticArrays: SVector
import MCMRSimulator.Plot: PlotPlane
import MCMRSimulator.Geometries.Internal: FixedSusceptibility, susceptibility_off_resonance
import MCMRSimulator.Geometries: fix_susceptibility


plot_off_resonance!(scene, plot_plane::PlotPlane, geometry) = plot_off_resonance!(scene, plot_plane, fix_susceptibility(geometry))
function plot_off_resonance!(scene, plot_plane::PlotPlane, susc::FixedSusceptibility; ngrid=400, colormap=:viridis)
    dims = -0.5:(1/ngrid):0.5
    xx_1d = dims * plot_plane.sizex
    yy_1d = dims * plot_plane.sizey
    pos_plane = broadcast(
        (x, y) -> SVector{3}([x, y, 0.]),
        reshape(xx_1d, length(xx_1d), 1),
        reshape(yy_1d, 1, length(yy_1d)),
    )
    pos_orig = inv(plot_plane.transformation).(pos_plane)
    field = map(p->susceptibility_off_resonance(susc, p), pos_orig)
    Makie.image!(scene, xx_1d, yy_1d, field, colormap=colormap)
end

end