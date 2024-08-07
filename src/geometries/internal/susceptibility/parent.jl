module Parent

import MRIBuilder: B0
import StaticArrays: SVector, SMatrix
import LinearAlgebra: norm
import .....Constants: gyromagnetic_ratio
import ...HitGrids: HitGrid, get_objects
import ...Obstructions: FixedObstruction 
import ...BoundingBoxes: BoundingBox 
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient


struct ShiftedSusceptibility{N, O<:BaseSusceptibility{N}} <: FixedObstruction{N}
    position :: SVector{N, Float64}
    radius :: Float64
    base :: O
end

BoundingBox(s::ShiftedSusceptibility) = BoundingBox(s.position .- s.radius, s.position .+ s.radius)


"""
    ParentSusceptibility(base::BaseSusceptibility, rotation, repeats, lorentz_radius)

Group of [`BaseSusceptibility`](@ref) susceptibility sources.
"""
struct ParentSusceptibility{N, O<:BaseSusceptibility{N}, R<:Union{Nothing, SVector{N, Float64}}, K}
    grid :: HitGrid{N, ShiftedSusceptibility{N, O}}
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
    susceptibility_off_resonance(susceptibility, position[, inside])

Computes the susceptibility off-resonance caused by a single [`BaseSusceptibility`](@ref) 
or all the susceptibility sources in a simulation or [`FixedSusceptibility`](@ref) at the given position.

The field is computed in ppm. Knowledge of the scanner [`B0`](@ref) is needed to convert it into KHz.
"""
function susceptibility_off_resonance(parent::ParentSusceptibility, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing) 
    rotated = parent.rotation * position
    if isnothing(parent.half_repeats)
        return susceptibility_off_resonance_non_repeating(parent, rotated, inside)
    else
        return susceptibility_off_resonance_repeating(parent, rotated, inside)
    end
end

function susceptibility_off_resonance_non_repeating(parent::ParentSusceptibility{N}, position::SVector{N, Float64}, inside::Union{Nothing, Bool}) where {N}
    field = zero(Float64)
    b0_field = parent.rotation[:, 3]
    for (_, shifted) in get_objects(parent.grid, position)
        offset = position .- shifted.position
        dist = norm(offset)
        if isinf(parent.lorentz_radius) || dist - shifted.radius < parent.lorentz_radius
            field += single_susceptibility(shifted.base, offset, dist, inside, b0_field)
        end
    end
    return field
end

function susceptibility_off_resonance_repeating(parent::ParentSusceptibility{N}, position::SVector{N, Float64}, inside::Union{Nothing, Bool}) where {N}
    field = zero(Float64)

    normed = @. mod(position + parent.half_repeats, 2 * parent.half_repeats) - parent.half_repeats
    b0_field = parent.rotation[:, 3]

    for (use_pos, shifted) in get_objects(parent.grid, normed)
        offset = use_pos .- shifted.center
        if N == 1
            dist = abs(offset[1])
        elseif N == 2
            dist = sqrt(offset[1] * offset[1] + offset[2] * offset[2])
        else
            dist = sqrt(offset[1] * offset[1] + offset[2] * offset[2] + offset[3] * offset[3])
        end
        if dist < (parent.lorentz_radius + shifted.radius)
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