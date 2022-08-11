"""
Supertype of any obstruction to the free diffusion of water (e.g., [`Wall`](@ref), [`Mesh`](@ref), [`Cylinder`](@ref), [`Sphere`](@ref)).
"""
abstract type Obstruction end

"Sequence of [`Obstruction`](@ref) objects that the particles have to navigate."
const Obstructions{N, T} = SVector{N, T} where {T <: Obstruction}

"""
    isinside(position, obstructions)
    isinside(spin, obstructions)
    isinside(snapshot, obstructions)

Test whether the particles are inside any of the [`Obstrunctions`](@ref) (or in the [`BoundingBox`](@ref)).
"""
isinside(pos::Vector, o) = isinside(PosVector(pos), o)
isinside(spin::Union{Spin, MultiSpin}, o) = isinside(position(spin), o)
isinside(snapshot::Union{Snapshot, MultiSnapshot}, o) = map(s -> isinside(s, o), snapshot.spins)
isinside(pos::PosVector, obstructions::Obstructions) = any(o -> isinside(pos, o), obstructions)


"""
    project(position, repeat::Repeated)
    project(position, transform::Transformed)

Computes the position in the space of the obstructions wrapped by the [`Repeated`](@ref) or [`Transformed`](@ref).

    project(position, grid::GridShape)

Computes the voxel index for the position on the [`GridShape`](@ref). 
This will return a result even if the point is outside of the grid. Use [`isinside`](@ref)(position, grid) to check that.
"""
function project end


include("bounding_box.jl")
include("diffuse.jl")
include("grid.jl")
include("repeat.jl")
include("transform.jl")
include("wall.jl")
include("sphere.jl")
include("cylinder.jl")
include("mesh.jl")