module OffResonance
using Makie
import StaticArrays: SVector
import MCMRSimulator.Plot: Plot, PlotPlane, plot_off_resonance!
import MCMRSimulator.Geometries.Internal: FixedSusceptibility, susceptibility_off_resonance
import MCMRSimulator.Geometries: fix_susceptibility
import ..Geometries: GeometryLike


Plot.plot_off_resonance!(scene, plot_plane::PlotPlane, geometry) = plot_off_resonance!(scene, plot_plane, fix_susceptibility(geometry))
function Plot.plot_off_resonance!(scene, plot_plane::PlotPlane, susc::FixedSusceptibility; ngrid=400, colormap=:viridis)
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

function Plot.plot_off_resonance(plot_plane::PlotPlane, geometry; figure=Dict{Symbol, Any}(), axis=Dict{Symbol, Any}(), kwargs...)
    f = Figure(; figure...)
    axis[:xgridvisible] = pop!(axis, :xgridvisible, false)
    axis[:ygridvisible] = pop!(axis, :ygridvisible, false)
    ax = Axis(f[1, 1]; axis...)
    plot_off_resonance!(ax, plot_plane, geometry; kwargs...)
    if length(ax.scene.plots) == 0
        return Makie.FigureAxis(f, ax)
    else
        return Makie.FigureAxisPlot(f, ax, ax.scene.plots[end])
    end
end

end