@Makie.recipe(SnapshotPlot, snap) do scene
    Makie.Theme(
    )
end


"""
    plot(snapshot)
    plot!(snapshot)

Plots the spin positions in the [`Snapshot`](@ref) in 3D color coded by the spin's orientation (see [`color`](@ref)).
"""
function Makie.plot!(sp::SnapshotPlot)
    snap = sp[1]
    colors = @lift color.($snap.spins)
    pos = @lift [Makie.Point3f(s.position) for s in $snap.spins]
    Makie.meshscatter!(sp, pos, color=colors)
    sp
end

Makie.plottype(::Union{MultiSnapshot, Snapshot}) = SnapshotPlot


@Makie.recipe(SnapshotPlanarPlot, snap, plane) do scene
    Makie.Attributes(
        markersize=0.1,
        marker=:circle,
        markerspace=:pixel,
        vector=0.1,
        arrowsize=0.1,
    )
end

"""
    plot(snapshot, plot_plane)
    plot!(snapshot, plot_plane)

Plots the spins in the [`Snapshot`](@ref) projected onto given [`PlotPlane`](@ref).
Each spin is represented by an arrow showing the transverse component of the spin.
The background color shows the average spin in the neighbourhood color coded as described in [`color`](@ref).
"""
function Makie.plot!(sp::SnapshotPlanarPlot)
    snap = sp[1]
    planar = sp[2]

    projection = @lift project_on_grid($planar, $snap)
    xs = @lift $projection[1]
    ys = @lift $projection[2]
    c = @lift color.($projection[3])
    Makie.image!(sp, xs, ys, c)

    pos = @lift [Makie.Point2f(s.position[1:2]) for s in transform($planar, $snap).spins]
    vl = sp[:vector]
    directions = @lift [Makie.Point2f(vector(s.orientation)[1:2] .* $vl) for s in $snap.spins]
    kwargs = Dict([sym => sp[sym] for sym in [:arrowsize]] )
    Makie.arrows!(sp, pos, directions; color=:black, kwargs...)
    sp
end

Makie.plottype(::Union{MultiSnapshot, Snapshot}, ::PlotPlane) = SnapshotPlanarPlot