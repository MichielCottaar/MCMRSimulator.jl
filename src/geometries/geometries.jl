"""
Defines the microstructural environment with which the spins interact.
"""
module Geometries
include("internal/internal.jl")
include("user/user.jl")
import .User: 
    fix, fix_susceptibility,
    ObstructionGroup, IndexedObstruction,
    Wall, Walls, walls,
    Cylinder, Cylinders, cylinders,
    Sphere, Spheres, spheres,
    Annulus, Annuli, annuli,
    Triangle, Mesh, mesh,
    load_mesh, random_positions_radii, split_mesh,
    write_geometry, read_geometry
end