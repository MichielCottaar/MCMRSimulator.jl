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
    plane = ss[1]
    snapshots = ss[2]
    nspins = @lift isnothing($(ss[:nspins])) ? length($snapshots[1]) : maximum(($(ss[:nspins]), length($snapshots[1])))
    for index in 1:nspins[]
        positions_3d = @lift map(s -> position(s[index]), $snapshots)
        projected = @lift project_trajectory($plane, $positions_3d)
        positions = @lift Makie.Point2f.($projected[1])
        colors_main = @lift map(s -> color(get_sequence(s[index], $(ss[:sequence]))), $snapshots)
        colors = @lift [isfinite(index) ? $colors_main[Int(round(index))] : Colors.HSV() for index in $projected[2]]
        Makie.lines!(ss, positions; color=colors)
    end
    ss
end

Makie.plottype(::PlotPlane, ::Vector{<:Snapshot}) = Plot_Trajectory2D