using GLMakie
using MRSimulator

## Create an interactive sequence plot
fig = Figure()
sg = SliderGrid(fig[2, 1],
	(label="TE", range=1:0.1:50, format="{:.1f}μs", startvalue=20),
)
TE = sg.sliders[1].value
sequence = @lift Sequence([RFPulse(flip_angle=90), RFPulse(flip_angle=180, time=$TE/2.)], 1.5 * $TE)
ax, _ = plot(fig[1, 1], sequence)
ylims!(ax, 0., nothing)
on(sequence) do s
	xlims!(ax, nothing, s.TR)
end


## Plot a single snapshot
snap = Snapshot([
	Spin(transverse=0.3), 
	Spin(position=[0, 1, 2.], transverse=1.),
	Spin(position=[2., 1, 2.], transverse=0.)], 0.3)
plot(snap)