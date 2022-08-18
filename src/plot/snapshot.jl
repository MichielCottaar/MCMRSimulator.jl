@Makie.recipe(Scatter_Snapshot, snap) do scene
    Makie.Theme(
    )
end


"""
    plot(snapshot)
    plot!(snapshot)
    scatter_snapshot(snapshot)
    scatter_snapshot!(snapshot)

Plots the spin positions in the [`Snapshot`](@ref) in 3D color coded by the spin's orientation (see [`color`](@ref)).
"""
function scatter_snapshot end

function Makie.plot!(ss::Scatter_Snapshot)
    snap = ss[1]
    colors = @lift color.($snap)
    pos = @lift Makie.Point3f.(mr.position.($snap))
    Makie.meshscatter!(ss, pos, color=colors)
    ss
end

Makie.plottype(::Snapshot) = Scatter_Snapshot


@Makie.recipe(Dyad_Snapshot, plane, snap) do scene
    Makie.Attributes(
        dyadlength=0.1,
        arrowsize=0.1,
        color=:black,
        sequence=1
    )
end

"""
    dyad_snapshot(plot_plane, snapshot; dyadlength=0.1, arrowsize=0.1, color=:black, sequence=1)
    dyad_snapshot(plot_plane, snapshot; dyadlength=0.1, arrowsize=0.1, color=:black, sequence=1)

Plots the spins in the [`Snapshot`](@ref) projected onto given [`PlotPlane`](@ref).
Each spin is represented by an arrow showing the transverse component of the spin.
"""
function Makie.plot!(sp::Dyad_Snapshot)
    planar = sp[1]
    snap = sp[2]

    pos = @lift [Makie.Point2f(s.position[1:2]) for s in project($planar, $snap)]
    vl = sp[:dyadlength]
    sequence = sp[:sequence]
    directions = @lift [Makie.Point2f(orientation(get_sequence(s, $sequence))[1:2] .* $vl) for s in $snap]
    kwargs = Dict([sym => sp[sym] for sym in [:arrowsize, :color]] )
    Makie.arrows!(sp, pos, directions; kwargs...)
    sp
end


@Makie.recipe(Image_Snapshot, plane, snap) do scene
    Makie.Attributes(
        sequence=1
    )
end

"""
    image_snapshot(plot_plane, snapshot; vectorlength=0.1, arrowsize=0.1, color=:black, sequence=1)
    image_snapshot!(plot_plane, snapshot; vectorlength=0.1, arrowsize=0.1, color=:black, sequence=1)

Plots the spins in the [`Snapshot`](@ref) projected onto given [`PlotPlane`](@ref).
The average spin orientation across the plot plane is plotted using the colour coding from [`color`](@ref).
"""
function image_snapshot end

function Makie.plot!(sp::Image_Snapshot)
    planar = sp[1]
    snap = sp[2]

    projection = @lift project_on_grid($planar, $snap)
    xs = @lift $projection[1]
    ys = @lift $projection[2]
    c = @lift color.($projection[3])
    Makie.heatmap!(sp, xs, ys, c; interpolate=true)
    sp
end


@Makie.recipe(Plot_Snapshot, plane, snap) do scene
    Makie.Attributes(
        ndyads=10,
        dyadlength=0.1,
        arrowsize=0.1,
        color=:black,
        sequence=1
    )
end

"""
    plot(plot_plane, snapshot; kwargs...)
    plot!(plot_plane, snapshot; kwargs...)
    plot_snapshot(plot_plane, snapshot; kwargs...)
    plot_snapshot!(plot_plane, snapshot; kwargs...)

Plots the spins in the [`Snapshot`](@ref) projected onto given [`PlotPlane`](@ref).
Each spin is represented by an arrow showing the transverse component of the spin (for a total of the first `ndyads` spins).
The average spin orientation across the plot plane is plotted using the colour coding from [`color`](@ref).

This combines the plotting from [`image_snapshot`](@ref) and [`dyad_snapshot`](@ref).
"""
function plot_snapshot end

function Makie.plot!(sp::Plot_Snapshot)
    planar = sp[1]
    snap = sp[2]

    kwargs = Dict([sym => sp[sym] for sym in [:sequence]] )
    image_snapshot!(sp, planar, snap; kwargs...)
    snap_dyads = @lift $snap[1:$(sp[:ndyads])]
    kwargs = Dict([sym => sp[sym] for sym in [:arrowsize, :color, :dyadlength, :sequence]] )
    dyad_snapshot!(sp, planar, snap_dyads; kwargs...)
    sp
end


Makie.plottype(::PlotPlane, ::Snapshot) = Plot_Snapshot