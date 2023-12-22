module Trajectory
using Makie
import Colors
import MCMRSimulator.Plot: PlotPlane, project_trajectory, Plot, plot_trajectory!
import ..Utils: color
import MCMRSimulator.Spins: Snapshot, position, get_sequence


function Plot.plot_trajectory(args...; figure=Dict{Symbol, Any}(), axis=Dict{Symbol, Any}(), kwargs...)
    f = Figure(; figure...)
    axis[:xgridvisible] = pop!(axis, :xgridvisible, false)
    axis[:ygridvisible] = pop!(axis, :ygridvisible, false)
    ax_type = length(args) == 1 ? Axis3 : Axis
    ax = ax_type(f[1, 1]; axis...)
    plot_trajectory!(ax, args...; kwargs...)
    if length(ax.scene.plots) == 0
        return Makie.FigureAxis(f, ax)
    else
        return Makie.FigureAxisPlot(f, ax, ax.scene.plots[end])
    end
end


function plot_trajectory!(scene, trajectory::Vector{<:Snapshot}; sequence=1, kwargs...)
    for index in 1:length(trajectory[1])
        positions = map(s -> Makie.Point3f(position(s[index])), trajectory)
        colors = map(s -> color(s[index]; sequence=sequence), trajectory)
        Makie.lines!(scene, positions; color=colors, kwargs...)
    end
end

function plot_trajectory!(scene, plot_plane::PlotPlane, trajectory::Vector{<:Snapshot}; sequence=1, kwargs...)
    function _get_colors(sequence_index :: Integer, spin_index :: Integer, snapshots :: Vector{Snapshot{N}}, times :: Vector{Float64}) where {N}
        if N == 0 || sequence_index == 0
            return [Colors.HSV() for _ in times]
        else
            colors_main = map(s -> color(get_sequence(s[spin_index], sequence_index)), snapshots)
            return [isfinite(time) ? colors_main[Int(round(time))] : Colors.HSV() for time in times]
        end
    end

    for index in 1:length(trajectory[1])
        positions_3d = map(s -> position(s[index]), trajectory)
        projected = project_trajectory(plot_plane, positions_3d)
        positions = Makie.Point2f.(projected[1])
        colors = _get_colors(sequence, index, trajectory, projected[2])
        Makie.lines!(scene, positions; color=colors)
    end
end

end