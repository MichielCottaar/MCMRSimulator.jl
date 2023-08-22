"""
User interface for defining geometry.
"""
module User
include("obstructions/obstructions.jl")
include("fix.jl")
include("fix_susceptibility.jl")
include("load_mesh.jl")
include("random_distribution.jl")
import .Fix: fix
import .FixSusceptibility: fix_susceptibility
import .Obstructions: 
    ObstructionGroup, IndexedObstruction,
    Wall, Walls, walls,
    Cylinder, Cylinders, cylinders,
    Sphere, Spheres, spheres,
    Annulus, Annuli, annuli,
    Triangle, Mesh, mesh
import .LoadMesh: load_mesh
import .RandomDistribution: random_positions_radii
end