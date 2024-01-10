module Snapshots
using Makie
import MCMRSimulator.Spins: Snapshot, get_sequence, position, orientation
import ..Utils: Utils
import MCMRSimulator.Plot: PlotPlane, project, project_on_grid, Plot_Snapshot


function Makie.plot!(scene::Plot_Snapshot)
    lift(scene[:kind]) do kind
        func = Dict(
            :scatter => scatter_snapshot!,
            :dyad => dyad_snapshot!,
            :image => image_snapshot!
        )[kind]
        func(scene)
    end
end

Makie.plottype(::Snapshot) = Plot_Snapshot
Makie.plottype(::PlotPlane, ::Snapshot) = Plot_Snapshot

Makie.convert_arguments(::Plot_Snapshot, pp::PlotPlane, snapshot::Snapshot) = (pp, snapshot)


function generic_kwargs(scene::Plot_Snapshot)
    return Dict([key => scene[key] for key in [
        :visible, :overdraw, :transparency, :fxaa, :inspectable, :depth_shift, :model, :space,
        :interpolate
    ]])
end

# 3-dimensional plotting
function scatter_snapshot!(scene::Plot_Snapshot{<:Tuple{<:Snapshot}})
    @extract scene (snapshot, sequence, color)
    kwargs = generic_kwargs(scene)
    colors = @lift $color == Makie.automatic ? Utils.color.($snapshot; sequence=$sequence) : $color
    pos = @lift Makie.Point3f.(position.($snapshot))
    Makie.meshscatter!(scene, pos; color=colors, kwargs...)
end

function dyad_snapshot!(scene::Plot_Snapshot{<:Tuple{<:Snapshot}})
    @extract scene (snapshot, sequence, color, dyadlength)
    kwargs = generic_kwargs(scene)
    colors = @lift $color == Makie.automatic ? Utils.color.($snapshot; sequence=$sequence) : $color
    pos = @lift Makie.Point3f.(position.($snapshot))
    directions = @lift [Makie.Point3f(orientation(get_sequence(s, $sequence)) .* $dyadlength) for s in $snapshot]
    Makie.arrows!(scene, pos, directions; color=colors, kwargs...)
end

function image_snapshot!(scene::Plot_Snapshot{<:Tuple{<:Snapshot}})
    error("3D plotting is not supported for snapshot plotting with kind=:image. Please select a different `kind` (:scatter or :dyad) or provide a PlotPlane.")
end


# 2-dimensional plotting
function scatter_snapshot!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot}})
    plot_plane = scene[1]
    snapshot = scene[2]
    @extract scene (sequence, color)
    kwargs = generic_kwargs(scene)
    colors = @lift $color == Makie.automatic ? Utils.color.($snapshot; sequence=$sequence) : $color
    pos = @lift [Makie.Point2f(project($plot_plane, position(spin))[1:2]) for spin in $snapshot]
    Makie.scatter!(scene, pos; color=colors, kwargs...)
end

function dyad_snapshot!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot}})
    plot_plane = scene[1]
    snapshot = scene[2]
    @extract scene (sequence, color, dyadlength)
    kwargs = generic_kwargs(scene)
    colors = @lift $color == Makie.automatic ? Utils.color.($snapshot; sequence=$sequence) : $color
    pos = @lift [Makie.Point2f(project($plot_plane, position(spin))[1:2]) for spin in $snapshot]
    directions = @lift [Makie.Point2f(orientation(get_sequence(s, $sequence))[1:2] .* $dyadlength) for s in $snapshot]
    Makie.arrows!(scene, pos, directions; color=colors, kwargs...)
end

function image_snapshot!(scene::Plot_Snapshot{<:Tuple{PlotPlane, Snapshot}})
    plot_plane = scene[1]
    snapshot = scene[2]
    @extract scene (sequence, ngrid)
    kwargs = generic_kwargs(scene)
    projection = @lift project_on_grid($plot_plane, get_sequence($snapshot, $sequence), $ngrid)
    @lift Makie.heatmap!(scene, $projection[1], $projection[2], Utils.color.($projection[3]); kwargs...)
end


end