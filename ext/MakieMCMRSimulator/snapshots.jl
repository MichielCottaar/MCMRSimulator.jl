module Snapshots
using Makie
import MCMRSimulator.Spins: Snapshot, get_sequence, position, orientation
import ..Utils: Utils
import MCMRSimulator.Plot: PlotPlane, project, project_on_grid, plot_snapshot, plot_snapshot!


function scatter_plot_attributes()
    Makie.@DocumentedAttributes begin
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
        arrowsize=Makie.automatic
        arrowhead=Makie.automatic
        arrowtail=Makie.automatic
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
        "Whether to interpolate for `kind=:image`."
        interpolate=true
        "Number of grid elements to plot for `kind=:image`."
        ngrid=20
    end
end

@recipe Plot_Snapshot (plot_plane::Union{Nothing, PlotPlane}, snapshot::Snapshot) begin
    "Plot format to use (:scatter, :dyad, or :image)."
    kind = :scatter
    "Set the color."
    color = Makie.automatic
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
Makie.convert_arguments(::Type{Plot_Snapshot}, pp::PlotPlane, snapshot::Snapshot) = (pp, snapshot)


@recipe Plot_Snapshot_Kind (plot_plane::Union{Nothing, PlotPlane}, snapshot::Snapshot, kind::Val) begin
    "Set the color."
    color = Makie.automatic
    "Which magnetisation to plot if multiple sequences were simulated."
    sequence = 1

    scatter_plot_attributes()...
    dyad_plot_attributes()...
    image_plot_attributes()...
    Makie.mixin_generic_plot_attributes()...
    ssao = false
end

Makie.argument_names(::Type{<: Plot_Snapshot_Kind}, N) = (:plot_plane, :snapshot, :val_kind)

function set_color_from_magnetisation!(scene)
    map!(scene.attributes, [:color, :snapshot, :sequence], :final_color) do color, snapshot, sequence
        if color == Makie.automatic
            return Utils.color.(snapshot; sequence=sequence)
        else
            return fill(color, length(snapshot))
        end
    end
end

function set_3d_position!(scene)
    map!(scene.attributes, [:snapshot], :position) do snapshot
        Makie.Point3f.(position.(snapshot))
    end
end

# 3-dimensional plotting
function Makie.plot!(scene::Plot_Snapshot_Kind{<:Tuple{Nothing, Snapshot, Val{:scatter}}})
    set_color_from_magnetisation!(scene)
    set_3d_position!(scene)
    Makie.scatter!(scene, scene.attributes, scene.position; color=scene.final_color)
end

function Makie.plot!(scene::Plot_Snapshot_Kind{<:Tuple{Nothing, Snapshot, Val{:dyad}}})
    set_color_from_magnetisation!(scene)
    set_3d_position!(scene)
    map!(scene.attributes, [:snapshot, :sequence], :directions) do snapshot, sequence
        [Makie.Point3f(orientation(get_sequence(s, sequence))) for s in snapshot]
    end
    Makie.arrows3d!(scene, scene.attributes, scene.position, scene.directions; color=scene.final_color)
end

function Makie.plot!(::Plot_Snapshot_Kind{<:Tuple{Nothing, <:Snapshot, Val{:image}}})
    error("3D plotting is not supported for snapshot plotting with kind=:image. Please select a different `kind` (:scatter or :dyad) or provide a PlotPlane.")
end


function set_2d_position!(scene)
    map!(scene.attributes, [:plot_plane, :snapshot], :position) do plot_plane, snapshot
        Makie.Point3f.(position.(snapshot))
        [Makie.Point2f(project(plot_plane, position(spin))[1:2]) for spin in snapshot]
    end
end

# 2-dimensional plotting
function Makie.plot!(scene::Plot_Snapshot_Kind{<:Tuple{PlotPlane, <:Snapshot, Val{:scatter}}})
    set_color_from_magnetisation!(scene)
    set_2d_position!(scene)
    Makie.scatter!(scene, scene.attributes, scene.position; color=scene.final_color)
end

function Makie.plot!(scene::Plot_Snapshot_Kind{<:Tuple{PlotPlane, <:Snapshot, Val{:dyad}}})
    set_color_from_magnetisation!(scene)
    set_2d_position!(scene)
    map!(scene.attributes, [:sequence, :snapshot], :directions) do sequence, snapshot
        # Transverse magnetisation is returned irrespective of plot plane orientation
        [Makie.Point2f(orientation(get_sequence(s, sequence))[1:2]) for s in snapshot]
    end
    Makie.arrows2d!(scene, scene.attributes, scene.position, scene.directions; color=scene.final_color)
end

function Makie.plot!(scene::Plot_Snapshot_Kind{<:Tuple{PlotPlane, <:Snapshot, Val{:image}}})
    Makie.register_computation!(scene.attributes, [:sequence, :plot_plane, :snapshot, :ngrid], [:x, :y, :matrix]) do inputs, changed, cached
        x, y, orientations = project_on_grid(inputs.plot_plane, get_sequence(inputs.snapshot, inputs.sequence), inputs.ngrid)
        return x, y, Utils.color.(orientations)
    end
    Makie.heatmap!(scene, scene.attributes, scene.x, scene.y, scene.matrix)
end


end
