"""
    plot_off_resonance(plot_plane, geometry)
    plot_off_resonance(plot_plane, geometry)

Plots the off-resonance of `geometry` in the [`PlotPlane`](@ref).
"""
@Makie.recipe(Plot_Off_Resonance, plot_plane, geometry) do scene
    Makie.Theme(
        colormap=:viridis
    )
end

function Makie.plot!(por::Plot_Off_Resonance)
    plot_plane = por[1]
    geometry = por[2]

    dims = @lift -0.5:(1/$plot_plane.ngrid):0.5
    xx_1d = @lift $dims * $plot_plane.sizex
    yy_1d = @lift $dims * $plot_plane.sizey
    pos_plane = @lift broadcast(
        (x, y) -> SVector{3}([x, y, 0.]),
        reshape($xx_1d, length($xx_1d), 1),
        reshape($yy_1d, 1, length($yy_1d)),
    )
    pos_orig = @lift inv($plot_plane.transformation).($pos_plane)
    field = @lift map(p->off_resonance($geometry, p), $pos_orig)
    Makie.image!(por, xx_1d, yy_1d, field, colormap=por[:colormap])
    por
end
