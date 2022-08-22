"""
Supertype of any obstruction to the free diffusion of water (e.g., [`Wall`](@ref), [`Mesh`](@ref), [`Cylinder`](@ref), [`Sphere`](@ref)).
"""
abstract type Obstruction end

"Sequence of [`Obstruction`](@ref) objects that the particles have to navigate."
const Obstructions{N, T} = SVector{N, T} where {T <: Obstruction}

"""
    isinside(position, obstructions/bounding_box)
    isinside(spin, obstructions/bounding_box)
    isinside(snapshot, obstructions/bounding_box)

Test whether the particles are inside any of the [`Obstrunctions`](@ref) (or in the [`BoundingBox`](@ref)).
"""
isinside(pos::Vector, o) = isinside(PosVector(pos), o)
isinside(spin::Spin, o) = isinside(position(spin), o)
isinside(snapshot::Snapshot, o) = map(s -> isinside(s, o), snapshot.spins)
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


"""
    off_resonance(obstructions, position[, b0_field])

Computes the off-resonance field at the position due to the obstructions with the magnetic field orientation from `b0_field`.
"""
function off_resonance(obstructions::Obstructions, position::PosVector, b0_field=PosVector([0, 0, 1])::PosVector)
    sum(o->off_resonance(o, position, b0_field), obstructions)
end


include("bounding_box.jl")
include("diffuse.jl")
include("grid.jl")
include("repeat.jl")
include("transform.jl")
include("wall.jl")
include("sphere.jl")
include("cylinder.jl")
include("mesh.jl")