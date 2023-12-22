module Plot
using Makie


include("utils.jl")
include("plot_planes.jl")
include("sequences.jl")
include("geometries.jl")
include("off_resonance.jl")
include("snapshots.jl")
include("trajectory.jl")
include("movie.jl")

import .PlotPlanes: PlotPlane
import .Geometries: plot_geometry, plot_geometry!, plot_geometry3d, plot_geometry3d!
import .Snapshots: plot_snapshot, image_snapshot, dyad_snapshot
import .Snapshots: plot_snapshot!, image_snapshot!, dyad_snapshot!
import .Trajectory: plot_trajectory2d, plot_trajectory3d
import .Trajectory: plot_trajectory2d!, plot_trajectory3d!
import .OffResonance: plot_off_resonance, plot_off_resonance!
import .Movie: simulator_movie
end