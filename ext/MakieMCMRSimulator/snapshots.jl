module Snapshots
using Makie
import MCMRSimulator.Spins: Snapshot, get_sequence, position, orientation
import ..Utils: color
import MCMRSimulator.Plot: Plot, PlotPlane, project, project_on_grid, plot_snapshot!

function Plot.plot_snapshot(args...; figure=Dict{Symbol, Any}(), axis=Dict{Symbol, Any}(), kwargs...)
    f = Figure(; figure...)
    axis[:xgridvisible] = pop!(axis, :xgridvisible, false)
    axis[:ygridvisible] = pop!(axis, :ygridvisible, false)
    ax_type = length(args) == 1 ? Axis3 : Axis
    ax = ax_type(f[1, 1]; axis...)
    plot_snapshot!(ax, args...; kwargs...)
    if length(ax.scene.plots) == 0
        return Makie.FigureAxis(f, ax)
    else
        return Makie.FigureAxisPlot(f, ax, ax.scene.plots[end])
    end
end

function Plot.plot_snapshot!(args...; kind=:scatter, kwargs...)
    func = Dict(
        :scatter => scatter_snapshot!,
        :dyad => dyad_snapshot!,
        :image => image_snapshot!,
    )[kind]
    func(args...; kwargs...)
end

# 3-dimensional plotting
function scatter_snapshot!(scene, snapshot::Snapshot; sequence=1, kwargs...)
    colors = color.(snap; sequence=sequence)
    pos = Makie.Point3f.(position.(snap))
    Makie.meshscatter!(scene, pos; color=colors, kwargs...)
end

function dyad_snapshot!(scene, snapshot::Snapshot; sequence=1, dyad_length=0.1, kwargs...)
    pos = Makie.Point3f.(position.(snap))
    directions = [Makie.Point3f(orientation(get_sequence(s, sequence)) .* dyad_length) for s in snap]
    Makie.arrows!(scene, pos, directions; kwargs...)
end

function image_snapshot!(scene, snapshot::Snapshot; kwargs...)
    error("3D plotting is not supported for snapshot plotting with kind=:image. Please select a different `kind` (:scatter or :dyad) or provide a PlotPlane.")
end


# 2-dimensional plotting
function scatter_snapshot!(scene, plot_plane::PlotPlane, snapshot::Snapshot; sequence=1, kwargs...)
    colors = color.(snap; sequence=sequence)
    pos = [Makie.Point2f(project(plot_plane, position(spin))[1:2]) for spin in snapshot]
    Makie.meshscatter!(scene, pos; color=colors, kwargs...)
end

function dyad_snapshot!(scene, plot_plane::PlotPlane, snapshot::Snapshot; sequence=1, dyad_length=0.1, kwargs...)
    pos = [Makie.Point2f(project(plot_plane, position(spin))[1:2]) for spin in snapshot]
    directions = [Makie.Point2f(orientation(get_sequence(s, sequence))[1:2] .* dyad_length) for s in snapshot]
    Makie.arrows!(scene, pos, directions; kwargs...)
end

function image_snapshot!(scene, plot_plane::PlotPlane, snapshot::Snapshot; sequence=1, ngrid=20, interpolate=true, kwargs...)
    projection = project_on_grid(planar, get_sequence(snapshot, sequence=sequence), ngrid)
    Makie.heatmap!(sp, projection[1], projection[2], color.(projection[3]); kwargs...)
end

end