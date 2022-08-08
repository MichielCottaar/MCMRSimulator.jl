"gyromagnetic ratio in 10^6 rad/ms/T."
const gyromagnetic_ratio = 0.26752218744  # (10^6 rad⋅ms^−1⋅T^−1)

"Float type used during the simulation (default: Float64). This has to be changed in the code."
const Float = Float64
const SA = SA_F64

"General definition used for length-3 vectors. Mostly used for positions."
const PosVector = SVector{3, Float}

"""
Supertype of any obstruction to the free diffusion of water (e.g., [`Wall`](@ref), mesh, [`Cylinder`](@ref), [`Sphere`](@ref)).
"""
abstract type Obstruction end

"Sequence of [`Obstruction`](@ref) objects that the particles have to navigate."
const Obstructions{N, T} = SVector{N, T} where {T <: Obstruction}

"Intermediate object used internally to represent a movement from one position to another"
struct Movement
    origin :: PosVector
    destination :: PosVector
    timestep :: Float
end
