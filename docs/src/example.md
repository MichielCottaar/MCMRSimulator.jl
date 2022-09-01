# Example of diffusion in myelinated axons
## Generating and plotting the myelinated axons
Plotting is enabled through one of the Makie backends.
Here we will use CairoMakie.
```@example 1
using MRSimulator
using CairoMakie
```

First we define the geometry.
For this we will use [`random_annuli`](@ref), which produces random annuli distributed through space with a target density of 0.75 repeating every 20 micrometer:
```@example 1
import Random; Random.seed!(1)  # hide
geometry = random_annuli(0.75, repeats=[20., 20.], g_ratio=0.6, myelin=true, rotation=:y)
```

Before we can visualise this geometry, we must define the plane in which to make this plot.
As the annuli have been defined to be oriented in the y-direction, we also place the normal
of this [`PlotPlane`](@ref) in the y-direction
```@example 1
pp = PlotPlane(:y, size=20.)  # Define the plane in which to make the plot
f = plot_off_resonance(pp, geometry)  # plot the off-resonance field produced by the annuli
plot!(pp, geometry)  # overlay the geometry
save("myelinated_annuli.png", f); # hide
```
![](myelinated_annuli.png)
