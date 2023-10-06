module Parent

import StaticArrays: SVector, SMatrix
import LinearAlgebra: norm
import .....Constants: gyromagnetic_ratio
import .....Scanners: B0
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient
import ...Gridify: Grid, get_indices

"""
    ParentSusceptibility(base::BaseSusceptibility, rotation, repeats, lorentz_radius)

Group of `L` [`BaseSusceptibility`](@ref) susceptibility sources.
"""
struct ParentSusceptibility{L, N, O<:BaseSusceptibility{N}, R<:Union{Nothing, SVector{N, Float64}}, K}
    base :: Vector{O}
    grid :: Grid{N}
    positions_radii :: Vector{Tuple{SVector{N, Float64}, Float64}}
    rotation :: SMatrix{N, 3, Float64, K}
    half_repeats :: R
    lorentz_radius :: Float64
    maximum_radius :: Float64
end

"""
Tuple containing all susceptibility sources
"""
const FixedSusceptibility{N} = NTuple{N, ParentSusceptibility}

"""
    susceptibility_off_resonance(susceptibility/geometry, position[, inside])

Computes the susceptibility off-resonance caused by a single [`BaseSusceptibility`](@ref) 
or all the susceptibility sources in a [`FixedGeometry`](@ref) at the given position.

The field is computed in ppm. Knowledge of the scanner [`B0`](@ref) is needed to convert it into KHz.
"""
function susceptibility_off_resonance(parent::ParentSusceptibility{L, N}, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing) where {L, N}
    rotated = parent.rotation * position
    if isnothing(parent.half_repeats)
        return susceptibility_off_resonance_non_repeating(parent, rotated, inside)
    else
        return susceptibility_off_resonance_repeating(parent, rotated, inside)
    end
end

function susceptibility_off_resonance_non_repeating(parent::ParentSusceptibility{L, N}, position::SVector{N, Float64}, inside::Union{Nothing, Bool}) where {L, N}
    field = zero(Float64)
    b0_field = parent.rotation[:, 3]
    for (index, _) in get_indices(parent.grid, position)
        (center, radius) = parent.positions_radii[index]
        offset = position - center
        dist = norm(offset)
        if isinf(parent.lorentz_radius) || dist - radius < parent.lorentz_radius
            field += single_susceptibility(parent.base[index], offset, dist, inside, b0_field)
        end
    end
    return field
end

function susceptibility_off_resonance_repeating(parent::ParentSusceptibility{L, N}, position::SVector{N, Float64}, inside::Union{Nothing, Bool}) where {L, N}
    field = zero(Float64)

    normed = @. mod(position + parent.half_repeats, 2 * parent.half_repeats) - parent.half_repeats
    b0_field = parent.rotation[:, 3]

    for (index, shift) in get_indices(parent.grid, normed)
        (center, radius) = parent.positions_radii[index]
        if iszero(shift)
            offset = normed .- center
        else
            offset = normed .- center .- parent.grid.shifts[shift]
        end
        if N == 1
            dist = abs(offset[1])
        elseif N == 2
            dist = sqrt(offset[1] * offset[1] + offset[2] * offset[2])
        else
            dist = sqrt(offset[1] * offset[1] + offset[2] * offset[2] + offset[3] * offset[3])
        end
        if dist < (parent.lorentz_radius + radius)
            field += single_susceptibility(parent.base[index], offset, dist, inside, b0_field)
        end
    end
    return field
end

function susceptibility_off_resonance(sources::FixedSusceptibility, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing)
    return sum([susceptibility_off_resonance(s, position, inside) for s in sources])
end

susceptibility_off_resonance(sources::FixedSusceptibility{0}, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing) = 0.

"""
    off_resonance_gradient(susceptibility, B0)

Maximum gradient of the off-resonance field in kHz/um due to the local susceptibility sources.
Internally, computed for each susceptibility sources using [`single_susceptibility_gradient`](@ref).
The maximum out of these is returned.
"""
function off_resonance_gradient(parent::ParentSusceptibility, B0)
    return maximum(single_susceptibility_gradient.(parent.base)) * B0 * gyromagnetic_ratio * 1e-6
end

function off_resonance_gradient(sources::FixedSusceptibility, B0)
    return maximum(off_resonance_gradient.(sources, B0))
end

off_resonance_gradient(sources::FixedSusceptibility{0}, B0) = 0.

end