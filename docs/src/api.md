```@meta
CurrentModule = MCMRSimulator
```

# [MCMRSimulator.jl public API](@id api)
This is the public API for [MCMRSimulator](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl).
For a more user-friendly introduction, click [here](@ref Introduction).

## Running simulations
```@docs
Simulation
readout
evolve
Subset
Snapshot
BoundingBox
Spin
MCMRSimulator.get_subset
MCMRSimulator.SpinOrientation
MCMRSimulator.SpinOrientationSum
MCMRSimulator.TimeStep
```

## Creating geometry
### Geometry types
```@docs
Annuli
Cylinders
Walls
Spheres
Mesh
BendyCylinder
```
### Geometry helper functions
```@docs
load_mesh
random_positions_radii
```

## Querying simulation output
```@docs
position
longitudinal
transverse
phase
orientation
isinside
stuck
stuck_to
get_time
```

## Plotting
```@docs
PlotPlane
plot_snapshot
plot_geometry
plot_trajectory
plot_off_resonance
simulator_movie
```

## Probing MCMR internals
```@docs
MCMRSimulator.gyromagnetic_ratio
MCMRSimulator.get_rotation
MCMRSimulator.get_readouts
MCMRSimulator.IndexedReadout
MCMRSimulator.susceptibility_off_resonance
MCMRSimulator.ObstructionGroup
MCMRSimulator.IndexedObstruction
MCMRSimulator.nvolumes
MCMRSimulator.fix
MCMRSimulator.fix_susceptibility
MCMRSimulator.FixedGeometry
MCMRSimulator.surface_relaxation
MCMRSimulator.surface_density
MCMRSimulator.dwell_time
MCMRSimulator.permeability
MCMRSimulator.GlobalProperties
MCMRSimulator.R1
MCMRSimulator.R2
MCMRSimulator.off_resonance
```