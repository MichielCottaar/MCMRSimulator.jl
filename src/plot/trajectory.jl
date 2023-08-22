module Trajectory
using Makie
import Colors
import ..PlotPlanes: PlotPlane, project_trajectory
import ..Utils: color
import ...Spins: Snapshot, position, get_sequence
"""
    plot(snapshots)
    plot!(snapshots)
    plot_trajectory3d(snapshots)
    plot_trajectory3d!(snapshots)

Plots trajectory of the spins in a sequence of [`Snapshot`](@ref) objects.
"""
@Makie.recipe(Plot_Trajectory3D, snapshots) do scene
    Makie.Theme(
        nspins=nothing,
        sequence=1,
    )
end

function Makie.plot!(ss::Plot_Trajectory3D)
    snapshots = ss[1]
    nspins = @lift isnothing($(ss[:nspins])) ? length($snapshots[1]) : maximum(($(ss[:nspins]), length($snapshots[1])))
    for index in 1:nspins[]
        positions = @lift map(s -> Makie.Point3f(position(s[index])), $snapshots)
        colors = @lift map(s -> color(get_sequence(s[index], $(ss[:sequence]))), $snapshots)
        Makie.lines!(ss, positions; color=colors)
    end
    ss
end

Makie.plottype(::Vector{<:Snapshot}) = Plot_Trajectory3D


"""
    plot(snapshots)
    plot!(snapshots)
    plot_trajectory3d(snapshots)
    plot_trajectory3d!(snapshots)

Plots trajectory of the spins in a sequence of [`Snapshot`](@ref) objects.
"""
@Makie.recipe(Plot_Trajectory2D, plane, snapshots) do scene
    Makie.Theme(
        nspins=nothing,
        sequence=1,
    )
end

function Makie.plot!(ss::Plot_Trajectory2D)
    function _get_colors(sequence_index :: Integer, spin_index :: Integer, snapshots :: Vector{Snapshot{N}}, times :: Vector{Float64}) where {N}
        if N == 0 || sequence_index == 0
            return [Colors.HSV() for _ in times]
        else
            colors_main = map(s -> color(get_sequence(s[spin_index], sequence_index)), snapshots)
            return [isfinite(time) ? colors_main[Int(round(time))] : Colors.HSV() for time in times]
        end
    end
    plane = ss[1]
    snapshots = ss[2]
    nspins = @lift isnothing($(ss[:nspins])) ? length($snapshots[1]) : maximum(($(ss[:nspins]), length($snapshots[1])))
    for index in 1:nspins[]
        positions_3d = @lift map(s -> position(s[index]), $snapshots)
        projected = @lift project_trajectory($plane, $positions_3d)
        positions = @lift Makie.Point2f.($projected[1])
        colors = @lift _get_colors($(ss[:sequence]), index, $snapshots, $projected[2])
        Makie.lines!(ss, positions; color=colors)
    end
    ss
end

Makie.plottype(::PlotPlane, ::Vector{<:Snapshot}) = Plot_Trajectory2D

end