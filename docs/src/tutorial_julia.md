# [Tutorial](@id tutorial_julia)
This tutorial will walk through an example of modelling the MRI signal evolution for a diffusion-weighted sequence.
The spins in this simulation will be constrained by regularly packed cylinders.
This tutorial will use the programatic Julia interface, which you can run in the Julia REPL or in a Jupyter notebook.
If you would prefer to use the command line interface, you can find a tutorial doing the same analysis [here](@ref tutotial_cli).
After [installation](@ref installation) we can load MCMRSimulator.jl using
```@example tutorial
using MCMRSimulator
using CairoMakie  # used for plotting; use GLMakie or WGLMakie for interactive plots
```

In general, running a simulation will consist of the following three steps:
- Defining the microstructure and one or more sequences by creating an appropriate [`Simulation`](@ref) object.
- Initialising [`Snapshot`](@ref) with one or more [`Spin`](@ref) objects.
- Simulating a random walk of the spins through the microstructure and the MR signal produced by those spins.
- Plotting the MR signal or storing it to disk.
We will look through each of these steps below.
## Defining the simulation
The first step is to define the environment through which the spins will evolve.
We will do so by creating an appropriate [`Simulation`](@ref) object.
This `Simulation` will contain information on the microstructure, how spins interact with that microstructure, and the enabled sequence(s).

These different steps are described in more detail in other sections of this documentation:
- [How to define the microstrutural geometry](@ref geometry)
- [Controlling spin behaviour](@ref properties)
- [Sequence generation](@ref sequence)

First we will define a geometry formed of regularly packed axons.
This is represented by a single cylinder with a radius of 1 micrometer that repeats itself every 2.5 micrometer (in both the x-, and y-direction).
```@example tutorial
import Random; Random.seed!(1) # hide
geometry = cylinders(radius=1., repeats=[2.5, 2.5])

f = plot(PlotPlane(size=5), geometry)
f
save("tutorial_geometry.png", f) # hide
nothing # hide
```
![](tutorial_geometry.png)

More complicated geometries can be generated as described [here](@ref geometry).

The next step is to define a sequence (see [here](@ref sequence) for more details). 
Here we will adopt a single diffusion-weighted MRI sequence.
```@example tutorial
sequence = dwi(bval=2., TR=300, TE=80, scanner=Siemens_Prisma)  # default gradient orientation in the x-direction
f = plot(sequence)
f
save("tutorial_sequence.png", f); # hide
nothing # hide
```
![](tutorial_sequence.png)

Once we have both a geometry and one or more sequences, we can put them together in a [`Simulation`](@ref) object:
```@example tutorial
simulation = Simulation(sequence, R2=0.012, R1=3e-3, diffusivity=2., off_resonance=0.1, geometry=geometry)
nothing # hide
```
By default there is no T1 or T2 relaxation and spins do not move.
Enabling spin relaxation and diffusion requires setting the appropriate parameters in the [`Simulation`](@ref).

## Initialising the simulation
The current state of the simulation at any time is given by a [`Snapshot`](@ref) object.
This is essentially a vector of [`Spin`](@ref) objects with a time stamp.
Each [`Spin`](@ref) represents a single diffusing particle.
Besides containing its current position, it also contains its contribution to the MR signal for each of the sequences in the simulation and whether it is stuck on any surfaces.

The recommended way to initialise is to call [`Snapshot`](@ref)`(<number of spins>, <simulation>, [bounding_box])`.
This will create randomly distributed spins within some [`BoundingBox`](@ref).
By default this bounding box is an isotropic voxel with a size of 1 mm centered on the origin.

After initialisation or after running the simulation, the [`Snapshot`](@ref) can be later filtered to, for example, only include intra-cellular or extracellular spins (see [`isinside`](@ref)) or only free or stuck spins (see [`stuck`](@ref)).

The simulation can also be initialised explicitly using a sequence of positions (i.e., length-3 vectors) with the initial spin positions. 
Note that such a simulation will start with all spins free and not necessarily randomly distributed, which means it might take some time to reach an equilibrium.

For each of these initialisations the initial magnetisation can be explicitly set using the [`transverse`](@ref), [`longitudinal`](@ref), and [`phase`](@ref) flags.
The default is for spins to start in equilibrium (i.e., transverse magnetisation of 0 and longitudinal magnetisation of 1).

Finally, one could start a simulation using a [`Snapshot`](@ref) from a previous simulation.

!!! note "Deterministic spins"
    Each [`Spin`](@ref) is assigned a random number state at creation, which will be used for its future evolution. This means that after creation of a spin or a [`Snapshot`](@ref) its future is fully determined. This ensures that when a spin is evolved through the same simulation multiple times, it sill follow the same path each time. This allows improved comparisons between simulations with the same geometry, but different sequences/physics. However, it can lead to confusing results (e.g., a simulation initialised with `fill(Spin(), 500)` will contain 500 spins all following the exact same path).
## Running the simulation
Running the simulation is done through 4 functions, for which examples are shown below:
- [`trajectory`](@ref): follow the full state evolution for a small number of spins
- [`evolve`](@ref): evolve a large number of spins for a specific time
- [`readout`](@ref): return the spin states at the time of the sequence readouts
- [`signal`](@ref): return the average signal at high temporal resolution
For each of these functions a [`Snapshot`](@ref) object can be passed on or a number. 
In case of a number an initial [`Snapshot`](@ref) is created on the fly with the appropriate number of spins.
### Illustrating trajectories
We will start by illustrating the 2D trajectory for two spins, one inside and one outside of the cylinder.
To plot the trajectory we first need to output the state of the all spins at a high temporal resolution,
which can be done using [`trajectory`](@ref):
```@example tutorial
# Simulate 2 spins with given starting positions for 3 ms
snapshots = trajectory([[0, 0, 0], [1, 1, 0]], simulation, 0:0.01:3.)

pp = PlotPlane(size=5.)
f = plot(pp, geometry)
plot_trajectory2d!(pp, snapshots)
f
save("tutorial_trajectory2D.png", f) # hide
nothing # hide
```
![](tutorial_trajectory2D.png)
In this plot the color at each timepoint encodes the spin orientation.
The brightness of the spin indicates the size of the transverse component with purely longitudinal spins being in black.
The color of the spin encodes the phase of the MR signal in the transverse plane.

The trajectories can also be plotted in 3D:
```@example tutorial
f = plot_trajectory3d(snapshots)
save("tutorial_trajectory3D.png", f) # hide
nothing # hide
```
![](tutorial_trajectory3D.png)

When simulating a large number of spins, storing the spin state every timestep would become very memory-intensive.
To get around this we will typically either store the average signal at a high temporal resolution
or the full state at a low resolution. Let's have a look of examples for both.
### Signal evolution
To store the total signal evolution we can use [`signal`](@ref).
At each timepoint this will return the total MR signal (for each sequence) as a [`SpinOrientation`](@ref) object.
From this one can estimate the [`transverse`](@ref) component, the [`longitudinal`](@ref) component, and the [`phase`](@ref).
Here we plot the transverse component of the signal evolution as an example
```@example tutorial
times = 0:0.1:sequence.TR
average_signals = signal(3000, simulation, times)  # simulate 3000 spins for a single repetition time
f = plot(sequence)
lines!(times, transverse.(average_signals)/3000.)
f
save("tutorial_transverse.png", f) # hide
nothing # hide
```
![](tutorial_transverse.png)

### Readout at specific times
We can return a [`Snapshot`](@ref) at any time simply by running:
```@example tutorial
snapshot = evolve(3000, simulation, 80.)
pp = PlotPlane(size=2.5)
f = plot(pp, snapshot)
plot!(pp, geometry)
f
save("tutorial_snapshot.png", f) # hide
nothing # hide
```
![](tutorial_snapshot.png)
The color encoding is the same as for the trajectory plot above.
The brightness encodes the size of the transverse component, while
the color encodes the phase of the MR signal in the transverse plane.
We can see that outside of the cylinder the signal contribution is significantly reduced.
The black arrows show the transverse spin for some random spins.

The snapshot returned by [`evolve`](@ref) can be used as a starting point for further simulations.
We can use this to plot the longitudinal signal at the first and third TR using:
```@example tutorial
first_TR_start = Snapshot(1000)
fifth_TR_start = evolve(first_TR_start, simulation, sequence.TR * 2)

f = plot(sequence)

times = 0:0.1:100.
for start in (first_TR_start, fifth_TR_start)
    simulated_signals = signal(start, simulation, times .+ start.time)
    lines!(times, longitudinal.(simulated_signals)/3000., cycle=[:color])
end
f
save("tutorial_longitudinal.png", f) # hide
nothing # hide
```
![](tutorial_longitudinal.png)

Sequences will usually contain one or more [`Readout`](@ref) objects to mark the readout times.
To get the [`Snapshot`](@ref) at these readouts during one repetition time, you can use [`readout`](@ref).
