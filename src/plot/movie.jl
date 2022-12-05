function simulator_movie(filename, simulator::Simulation{N}, times, repeats; resolution=(1600, 800), trajectory_init=30, signal_init=10000, framerate=50) where {N}
    if isa(trajectory_init, Integer)
        trajectory_init = [rand(3) .* repeats .- repeats ./ 2 for _ in 1:trajectory_init]
    end
    sig = signal(signal_init, simulator, times)
    trans = [transverse.([s[index] for s in sig]) ./ signal_init for index in 1:N]
    traj = trajectory(trajectory_init, simulator, times);
    pp = PlotPlane(sizex=repeats[1], sizey=repeats[2])

    index = Observable(1)
    time = @lift times[$index]
    fig = Figure(resolution=resolution)
    ax_walk = Axis(fig[1, 1], title=(@lift("time = $(round($time, digits=1)) ms")))
    for geometry in simulator.micro.geometry
        plot!(ax_walk, pp, geometry)
    end
    xlims!(ax_walk, -5, 5)
    ylims!(ax_walk, -5, 5)
    dyad_snapshot!(ax_walk, pp, @lift(traj[$index]), dyadlength=0.6, arrowsize=10.)

    ax_seq = Axis(fig[1, 2])

    function set_nan(arr, index)
        new = copy(arr)
        new[index+1:end] .= NaN
        new
    end

    plot!(ax_seq, simulator.sequences[1])
    xlims!(ax_seq, -2, maximum(times) + 2)
    for t in trans
        lines!(ax_seq, times, (@lift(set_nan(t, $index))))
    end


    record(fig, filename, 1:length(times);
            framerate=framerate) do i
        index[] = i
    end
end