"""
Defines the microstructural environment with which the spins interact.
"""
module Geometries
include("utils.jl")
include("bounding_box.jl")
include("internal/internal.jl")
include("user/user.jl")
include("fix/fix.jl")
import .Fix:
    fix, fix_susceptibility
import .Internal: FixedGeometry, Intersection, ObstructionIndex
import .BoundingBoxes: BoundingBox
import .User:
    ObstructionGroup, IndexedObstruction,
    Wall, Walls,
    Cylinder, Cylinders,
    Sphere, Spheres,
    Annulus, Annuli,
    Triangle, Mesh,
    Ring, BendyCylinder,
    load_mesh, SWCFile, SWCNode, read_swc, read_swc_raw,
    random_positions_radii, nvolumes,
    write_geometry, read_geometry_json, read_geometry
end
