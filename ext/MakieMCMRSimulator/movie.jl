module Movie
using Makie
import MCMRSimulator.Plot: PlotPlane
import ..Geometries: plot_geometry!
import ..Snapshots: plot_snapshot!
import MCMRSimulator.Spins: transverse
import MCMRSimulator.Simulations: Simulation
import MCMRSimulator.Evolve: readout

function simulator_movie(filename, simulator::Simulation{N}, times, repeats; resolution=(1600, 800), trajectory_init=30, signal_init=10000, framerate=50, plane_orientation=:z, kwargs...) where {N}
    if isa(trajectory_init, Integer)
        trajectory_init = [rand(3) .* repeats .- repeats ./ 2 for _ in 1:trajectory_init]
    end
    sig = readout(signal_init, simulator, times)
    if ndims(sig) == 1
        trans = [transverse.(sig) ./ signal_init]
    else
        trans = [transverse.(sig[index, :]) ./ signal_init for index in 1:N]
    end
    traj = readout(trajectory_init, simulator, times, return_snapshot=true);
    if ndims(traj) > 1
        traj = traj[1, :]
    end
    pp = PlotPlane(plane_orientation, sizex=repeats[1], sizey=repeats[2])

    fig = Figure(resolution=resolution)
    function single_frame!(index::Integer)
        time = times[index]
        ax_walk = Axis(fig[1, 1], title="time = $(round(time, digits=1)) ms")
        plot_geometry!(ax_walk, pp, simulator.geometry)
        xlims!(ax_walk, -repeats[1]/2, repeats[1]/2)
        ylims!(ax_walk, -repeats[2]/2, repeats[2]/2)
        plot_snapshot!(ax_walk, pp, traj[index]; kwargs...)

        ax_seq = Axis(fig[1, 2])

        function set_nan(arr, index)
            new = copy(arr)
            new[index+1:end] .= NaN
            new
        end

        plot!(ax_seq, simulator.sequences[1])
        for t in trans
            lines!(ax_seq, times, set_nan(t, index))
        end
        xlims!(ax_seq, -3, maximum(times) + 5)
    end


    record(fig, filename, 1:length(times);
            framerate=framerate) do i
        empty!(fig)
        single_frame!(i)
    end
end

end