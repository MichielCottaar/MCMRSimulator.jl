"Intermediate object used internally to represent a movement from one position to another"
struct Movement{N}
    origin :: SVector{N, Float}
    destination :: SVector{N, Float}
    timestep :: Float
end

function Movement(origin::AbstractArray, destination::AbstractArray, timestep::Real=1) 
    ndim = length(origin)
    Movement{ndim}(SVector{ndim, Float}(origin), SVector{ndim, Float}(destination), Float(timestep))
end