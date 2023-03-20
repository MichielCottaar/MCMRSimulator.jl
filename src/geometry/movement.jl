"Intermediate object used internally to represent a movement from one position to another"
struct Movement{N}
    origin :: SVector{N, Float}
    destination :: SVector{N, Float}
end

function Movement(origin::AbstractArray, destination::AbstractArray) 
    ndim = length(origin)
    Movement{ndim}(SVector{ndim, Float}(origin), SVector{ndim, Float}(destination))
end