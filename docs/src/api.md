```@meta
CurrentModule = MCMRSimulator
```

# [MCMRSimulator.jl API](@id api)
This is the API for [MCMRSimulator](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl).
For a more user-friendly introduction, click [here](@ref Introduction).

```@index
```

```@autodocs
Modules = [
    MCMRSimulator,
    MCMRSimulator.Constants,
    MCMRSimulator.Methods,
    MCMRSimulator.Scanners,
    MCMRSimulator.Properties,
    MCMRSimulator.Geometries,
    MCMRSimulator.Geometries.User,
    MCMRSimulator.Geometries.User.Obstructions,
    MCMRSimulator.Geometries.User.Obstructions.Fields,
    MCMRSimulator.Geometries.User.Obstructions.ObstructionTypes,
    MCMRSimulator.Geometries.User.Obstructions.ObstructionGroups,
    MCMRSimulator.Geometries.User.Fix,
    MCMRSimulator.Geometries.User.FixSusceptibility,
    MCMRSimulator.Geometries.User.LoadMesh,
    MCMRSimulator.Geometries.User.RandomDistribution,
    MCMRSimulator.Geometries.Internal,
    MCMRSimulator.Geometries.Internal.RayGridIntersection,
    MCMRSimulator.Geometries.Internal.BoundingBoxes,
    MCMRSimulator.Geometries.Internal.Obstructions,
    MCMRSimulator.Geometries.Internal.Obstructions.ObstructionIntersections,
    MCMRSimulator.Geometries.Internal.Obstructions.FixedObstructions,
    MCMRSimulator.Geometries.Internal.Obstructions.Walls,
    MCMRSimulator.Geometries.Internal.Obstructions.Rounds,
    MCMRSimulator.Geometries.Internal.Obstructions.Triangles,
    MCMRSimulator.Geometries.Internal.Obstructions.Shifts,
    MCMRSimulator.Geometries.Internal.Gridify,
    MCMRSimulator.Geometries.Internal.Intersections,
    MCMRSimulator.Geometries.Internal.Reflections,
    MCMRSimulator.Geometries.Internal.FixedObstructionGroups,
    MCMRSimulator.Geometries.Internal.Properties,
    MCMRSimulator.Geometries.Internal.Susceptibility,
    MCMRSimulator.Geometries.Internal.Susceptibility.Base,
    MCMRSimulator.Geometries.Internal.Susceptibility.Parent,
    MCMRSimulator.Geometries.Internal.Susceptibility.Cylinder,
    MCMRSimulator.Geometries.Internal.Susceptibility.Annulus,
    MCMRSimulator.Spins,
    MCMRSimulator.Sequences,
    MCMRSimulator.Sequences.Methods,
    MCMRSimulator.Sequences.Instants,
    MCMRSimulator.Sequences.Shapes,
    MCMRSimulator.Sequences.Gradients,
    MCMRSimulator.Sequences.RadioFrequency,
    MCMRSimulator.Sequences.Main,
    MCMRSimulator.Sequences.PulseQ,
    MCMRSimulator.SequenceBuilder,
    MCMRSimulator.SequenceBuilder.BuildingBlocks,
    MCMRSimulator.SequenceBuilder.DefineSequence,
    MCMRSimulator.SequenceBuilder.Diffusion,
    MCMRSimulator.SequenceBuilder.Sequences,
    MCMRSimulator.SequenceBuilder.Sequences.GradientEcho,
    MCMRSimulator.SequenceBuilder.Sequences.SpinEcho,
    MCMRSimulator.Timestep,
    MCMRSimulator.Relax,
    MCMRSimulator.Simulations,
    MCMRSimulator.Evolve,
    MCMRSimulator.Plot,
    MCMRSimulator.Plot.Utils,
    MCMRSimulator.Plot.PlotPlanes,
    MCMRSimulator.Plot.Sequences,
    MCMRSimulator.Plot.Geometries,
    MCMRSimulator.Plot.OffResonance,
    MCMRSimulator.Plot.Snapshots,
    MCMRSimulator.Plot.Trajectory,
    MCMRSimulator.Plot.Movie,
]
```