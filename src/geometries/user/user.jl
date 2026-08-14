"""
User interface for defining geometry.
"""
module User
include("obstructions/obstructions.jl")
include("split_mesh.jl")
include("load_mesh.jl")
include("load_swc.jl")
include("random_distribution.jl")
include("json.jl")
include("load_geometry.jl")
include("to_mesh.jl")
import .Obstructions: ObstructionGroup, IndexedObstruction, nvolumes,
    Wall, Walls,
    Cylinder, Cylinders,
    Sphere, Spheres,
    Annulus, Annuli,
    Triangle, Mesh,
    Ring, BendyCylinder
import .RandomDistribution: random_positions_radii
import .LoadMesh: load_mesh
import .LoadSWC: SWCFile, SWCNode, read_swc, read_swc_raw
import .LoadGeometry: read_geometry
import .JSON: write_geometry, read_geometry_json
end
