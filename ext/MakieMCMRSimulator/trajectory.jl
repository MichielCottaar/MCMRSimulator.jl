module Trajectory
using Makie
import Colors
import MCMRSimulator.Plot: PlotPlane, project_trajectory
import ..Utils: Utils
import MCMRSimulator.Spins: Snapshot, position, get_sequence

@recipe Plot_Trajectory (plot_plane::Union{Nothing, PlotPlane}, trajectory::AbstractVector{<:Snapshot}) begin
    "Which magnetisation to plot if multiple sequences were simulated."
    sequence=1
    color=Makie.automatic
    Makie.mixin_generic_plot_attributes()...
end

Makie.argument_names(::Type{<: Plot_Trajectory}, N) = (:plot_plane, :trajectory)

function Makie.plot!(scene::Plot_Trajectory{<:Tuple{Nothing, <:AbstractVector{<:Snapshot}}})
    Makie.register_computation!(scene.attributes, [:trajectory, :color, :sequence], [:positions, :colors]) do inputs, changed, cached
        positions = Makie.Point3f[]
        colors = Any[]
        nspins = length(inputs.trajectory[1])
        for index in 1:nspins
            append!(positions, map(s -> Makie.Point3f(position(s[index])), inputs.trajectory))
            if inputs.color == Makie.automatic
                append!(colors, map(s -> Utils.color(s[index]; sequence=inputs.sequence), inputs.trajectory))
            else
                append!(colors, fill(inputs.color, length(inputs.trajectory)))
            end
            push!(positions, Makie.Point3f(NaN, NaN, NaN))
            push!(colors, Colors.HSV())
        end
        return positions, colors
    end
    Makie.lines!(scene, scene.attributes, scene.positions; color=scene.colors)
end

Makie.plottype(::AbstractVector{<:Snapshot}) = Plot_Trajectory


function Makie.plot!(scene::Plot_Trajectory{<:Tuple{<:PlotPlane, <:AbstractVector{<:Snapshot}}})
    function _get_colors(sequence_index :: Integer, spin_index :: Integer, snapshots :: Vector{<:Snapshot{N}}, times :: Vector{Float64}) where {N}
        if N == 0 || sequence_index == 0
            return [Colors.HSV() for _ in times]
        else
            colors_main = map(s -> Utils.color(get_sequence(s[spin_index], sequence_index)), snapshots)
            return [isfinite(time) ? colors_main[Int(round(time))] : Colors.HSV() for time in times]
        end
    end
    Makie.register_computation!(scene.attributes, [:plot_plane, :trajectory, :color, :sequence], [:positions, :colors]) do inputs, changed, cached
        positions = Makie.Point2f[]
        colors = Any[]
        nspins = length(inputs.trajectory[1])
        for index in 1:nspins
            positions_3d = map(s -> position(s[index]), inputs.trajectory)
            projected = project_trajectory(inputs.plot_plane, positions_3d)
            append!(positions, Makie.Point2f.(projected[1]))
            if inputs.color == Makie.automatic
                append!(colors, _get_colors(inputs.sequence, index, inputs.trajectory, projected[2]))
            else
                append!(colors, fill(inputs.color, length(inputs.trajectory)))
            end
            push!(positions, Makie.Point2f(NaN, NaN))
            push!(colors, Colors.HSV())
        end
        return positions, colors
    end
    Makie.lines!(scene, scene.attributes, scene.positions; color=scene.colors)
end

Makie.convert_arguments(::Type{Plot_Trajectory}, trajectory::AbstractVector{<:Snapshot}) = (nothing, trajectory)

Makie.plottype(::PlotPlane, ::AbstractVector{<:Snapshot}) = Plot_Trajectory

end
