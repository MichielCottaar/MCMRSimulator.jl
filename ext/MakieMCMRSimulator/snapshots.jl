module Snapshots
using Makie
import MCMRSimulator.Spins: Snapshot, get_sequence, position, orientation
import ..Utils: Utils
import MCMRSimulator.Plot: PlotPlane, project, project_on_grid, plot_snapshot, plot_snapshot!


function scatter_plot_attributes()
    Makie.@DocumentedAttributes begin
        "Used in 2D and 3D scatter plot of snapshots"
        marker = @inherit marker
        markersize = @inherit markersize
        strokecolor = @inherit markerstrokecolor
        strokewidth = @inherit markerstrokewidth
        glowcolor = (:black, 0.0)
        glowwidth = 0.0
    end
end

function dyad_plot_attributes()
    Makie.@DocumentedAttributes begin
        "Used in 2D and 3D dyad plot of snapshots"
        arrowsize=automatic
        arrowhead=automatic
        arrowtail=automatic
        linestyle=nothing
        lengthscale=1.
        quality=32
        markerspace=:pixel
        diffuse=0.4
        specular=0.2
        shininess=32f0
    end
end

function image_plot_attributes()
    Makie.@DocumentedAttributes begin
        "Used in 2D image plot of snapshots"
        "Whether to interpolate for `kind=:image`."
        interpolate=true
        "Number of grid elements to plot for `kind=:image`."
        ngrid=20
    end
end

"""
    plot([plot_plane], snapshot; kind=:scatter, kwargs...)
    plot!([scene,] [plot_plane], snapshot; kind=:scatter, kwargs...)
    plot_snapshot([plot_plane], snapshot; kind=:scatter, kwargs...)
    plot_snapshot!([scene,] [plot_plane], snapshot; kind=:scatter, kwargs...)

Plots a [`Snapshot`](@ref) in a new plot. 

The spins are plotted in 2D projected onto the [`PlotPlane`](@ref) if one is provided.
Otherwise, the spins are plotted in 3D (does not work for `kind=:image`).

There are three kinds of snapshot plots available:
## Scatter plot
Default (or set using `kind=:scatter`). Each spin is plotted as a point with the colour set by the transverse magnetisation.
Additional keywords are passed on to `Makie.scatter`
(namely, `marker`, `markersize`, `strokecolor`, `strokewidth`, `glowcolor`, `glowwidth`).

## Dyad plot
Set using `kind=:dyad`. Each spin is plotted as a dyad. For a 2D dyad the orienation is set by the transverse magnetisation.
For a 3D dyad the full magnetisation is used to set the orienation.
Additional keywords are passed on to `Makie.arrows`
(namely, `arrowsize`, `arrowhead`, `arrowtail`, `linestyle`, `lengthscale`, `quality`, `markerspace`, `diffuse`, `specular`, `shininess`).

## Image
Set using `kind=:image`. The average magnetisation is plotted across the [`PlotPlane`](@ref). 
The colour in each pixel is set by the average transverse magnetisation of the local spins.
Set `interpolate=false` to disable interpolating the colors.
This plot will not work in 3D (i.e., a [`PlotPlane`](@ref) is required).

This function will only work if [`Makie`](https://makie.org) is installed and imported.
"""
@recipe Plot_Snapshot (plot_plane::Union{Nothing, PlotPlane}, snapshot::Snapshot) begin
    "Plot format to use (:scatter, :dyad, or :image)."
    kind = :scatter
    "Set the color."
    color = automatic
    "Which magnetisation to plot if multiple sequences were simulated."
    sequence = 1

    scatter_plot_attributes()...
    dyad_plot_attributes()...
    image_plot_attributes()...
    Makie.mixin_generic_plot_attributes()...
    ssao = false
end

Makie.argument_names(::Type{<: Plot_Snapshot}, N) = (:plot_plane, :snapshot)

function Makie.plot!(scene::Plot_Snapshot)
    map!(scene.attributes, [:kind], :val_kind) do kind
        @assert kind in (:scatter, :dyad, :image) "Kind should be one of `:scatter`, `:dyad`, or `:image` for plotting snapshots, not `$kind`"
        Val(kind)
    end
    plot_snapshot_kind!(scene, scene.attributes, scene.plot_plane, scene.snapshot, scene.val_kind)
end

Makie.plottype(::Snapshot) = Plot_Snapshot
Makie.plottype(::PlotPlane, ::Snapshot) = Plot_Snapshot

Makie.convert_arguments(::Type{Plot_Snapshot}, snap::Snapshot) = (nothing, snap)
Makie.convert_arguments(::Plot_Snapshot, pp::PlotPlane, snapshot::Snapshot) = (pp, snapshot)


@recipe Plot_Snapshot_Kind (plot_plane::Union{Nothing, PlotPlane}, snapshot::Snapshot, kind::Val) begin
    "Set the color."
    color = automatic
    "Which magnetisation to plot if multiple sequences were simulated."
    sequence = 1

    scatter_plot_attributes()...
    dyad_plot_attributes()...
    image_plot_attributes()...
    mixin_generic_plot_attributes()...
    ssao = false
end

Makie.argument_names(::Type{<: Plot_Snapshot_Kind}, N) = (:plot_plane, :snapshot, :val_kind)

function set_color_from_magnetisation!(scene)
    map!(scene.attributes, [:color, :snapshot, :sequence], :final_color) do color, snapshot, sequence
        if color == Makie.automatic
            return color
        else
            return Utils.color.(snapshot; sequence=sequence)
        end
    end
end

function set_3d_position!(scene)
    map!(scene.attributes, [:snapshot], :position) do snapshot
        Makie.Point3f.(position.(snapshot))
    end
end

# 3-dimensional plotting
function plot_snapshot_kind!(scene::Plot_Snapshot_Kind{<:Tuple{Nothing, Snapshot, Val{:scatter}}})
    set_color_from_magnetisation!(scene)
    set_3d_position!(scene)
    Makie.scatter!(scene, scene.attributes, scene.position; color=scene.final_color)
end

function plot_snapshot_kind!(scene::Plot_Snapshot{<:Tuple{Nothing, Snapshot, Val{:dyad}}})
    set_color_from_magnetisation!(scene)
    set_3d_position!(scene)
    map!(scene.attributes, [:snapshot, :sequence], :directions) do snapshot, sequence
        [Makie.Point3f(orientation(get_sequence(s, sequence))) for s in snapshot]
    end
    Makie.arrows!(scene, scene.attributes, scene.position, directions; color=scene.final_color)
end

function plot_snapshot_kind!(::Plot_Snapshot_Kind{<:Tuple{Nothing, Snapshot, Val{:image}}})
    error("3D plotting is not supported for snapshot plotting with kind=:image. Please select a different `kind` (:scatter or :dyad) or provide a PlotPlane.")
end


function set_2d_position!(scene)
    map!(scene.attributes, [:plot_plane, :snapshot], :position) do plot_plane, snapshot
        Makie.Point3f.(position.(snapshot))
        [Makie.Point2f(project(plot_plane, position(spin))[1:2]) for spin in snapshot]
    end
end

# 2-dimensional plotting
function plot_snapshot_kind!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot, Val{:scatter}}})
    set_color_from_magnetisation!(scene)
    set_2d_position!(scene)
    Makie.scatter!(scene, scene.attributes, scene.position; color=scene.final_color)
end

function plot_snapshot_kind!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot, Val{:dyad}}})
    set_color_from_magnetisation!(scene)
    set_2d_position!(scene)
    map!(scene.attributes, [:sequence, :snapshot], :directions) do sequence, snapshot
        # Transverse magnetisation is returned irrespective of plot plane orientation
        [Makie.Point2f(orientation(get_sequence(s, sequence))[1:2]) for s in snapshot]
    end
    Makie.arrows!(scene, scene.attributes, scene.position, scene.directions; color=scene.final_color)
end

function plot_snapshot_kind!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot, Val{:image}}})
    Makie.register_computation!(scene.attributes, [:sequence, :plot_plane, :snapshot, :ngrid], [:x, :y, :matrix]) do inputs, changed, cached
        return project_on_grid(inputs.plot_plane, get_sequence(inputs.snapshot, inputs.sequence), inputs.ngrid)
    end
    Makie.heatmap!(scene, scene.attributes, scene.x, scene.y, scene.matrix)
end


end