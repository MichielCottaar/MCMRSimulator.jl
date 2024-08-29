module Grid

import StaticArrays: SVector, SMatrix
import LinearAlgebra: norm
import .....Constants: gyromagnetic_ratio
import ...BoundingBoxes: BoundingBox, lower
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient


struct SusceptibilityGridElement{N}
    index :: Int32
    shift :: Int32
    position :: SVector{N, Float64}
    radius :: Float64
    susceptibility :: Float64
end


BoundingBox(element::SusceptibilityGridElement{N}) where {N} = BoundingBox{N}(
    element.position .- radius,
    element.position .+ radius,
)

"""
Precomputed off-resonance field from repeating distant sources.

In each voxel the pre-computed field is stored in `off_resonance`.
The sources for which the field contribution should still be added are stored in `indices`.
These are stored of vectors of [`SuscetibilityGridElement`](@ref), where each element
contains all the information to compute the dipole field approximation.
"""
abstract type SusceptibilityGrid{N, O, K} end

"""
Specialised version of [`SusceptibilityGrid`](@ref) for repeating geometries.
"""
struct SusceptibilityGridRepeat{N, O, K} <: SusceptibilityGrid{N, O, K}
    inv_resolution :: SVector{N, Float64}
    rotation :: SMatrix{N, 3, Float64, K}

    off_resonance :: Array{Float64, N}
    half_repeats :: SVector{N, Float64}

    indices :: Array{Vector{SusceptibilityGridElement{N}}, N}
    sources :: Vector{O}
    shifts :: Vector{SVector{N, Float64}}
    B0_field :: SVector{N, Float64}
end

"""
Specialised version of [`SusceptibilityGrid`](@ref) for non-repeating geometries.
"""
struct SusceptibilityGridNoRepeat{N, O, K} <: SusceptibilityGrid{N, O, K}
    inv_resolution :: SVector{N, Float64}
    rotation :: SMatrix{N, 3, Float64, K}

    bounding_box_off_resonance :: BoundingBox{N}
    off_resonance :: Array{Float64, N}

    bounding_box_indices :: BoundingBox{N}
    indices :: Array{Vector{SusceptibilityGridElement{N}}, N}
    sources :: Vector{O}
    shifts :: Vector{SVector{N, Float64}}
    B0_field :: SVector{N, Float64}
    total_susceptibility :: Float64
    center_susceptibility :: SVector{N, Float64}
end

"""
Tuple containing all susceptibility sources
"""
const FixedSusceptibility{N} = NTuple{N, SusceptibilityGrid}


norm_position(::SusceptibilityGridNoRepeat, position) = position
norm_position(grid::SusceptibilityGridRepeat, position) = @. mod(position + grid.half_repeats, 2 * grid.half_repeats) - grid.half_repeats

get_coordinates(grid::SusceptibilityGrid, position, lower) = @. Int(floor((position - lower) * grid.inv_resolution)) + 1
get_coordinates(grid::SusceptibilityGridNoRepeat, position) = (
    get_coordinates(grid, position, lower(grid.bounding_box_off_resonance)),
    get_coordinates(grid, position, lower(grid.bounding_box_indices)),
)
function get_coordinates(grid::SusceptibilityGridRepeat, position)
    c = get_coordinates(grid, position, -grid.half_repeats)
    return (c, c)
end

"""
    susceptibility_off_resonance(susceptibility_grid, position[, inside])

Computes the susceptibility off-resonance caused by a [`SusceptibilityGrid`](@ref) at given position.

The field is computed in ppm. Knowledge of the scanner [`B0`](@ref) is needed to convert it into KHz.
"""
function susceptibility_off_resonance(grid::SusceptibilityGrid, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing)
    normed = norm_position(grid, grid.rotation * position)

    coord_off_resonance, coord_indices = get_coordinates(grid, normed)
    if any(coord_off_resonance .< 1) || any(coord_off_resonance .> size(grid.off_resonance))
        @assert grid isa SusceptibilityGridNoRepeat
        offset = normed - grid.center_susceptibility
        return dipole_approximation(grid.total_susceptibility, offset, norm(offset), grid.B0_field)
    end

    field = grid.off_resonance[coord_off_resonance...]
    if any(coord_indices .< 1) || any(coord_indices .> size(grid.indices))
        @assert grid isa SusceptibilityGridNoRepeat
        return field
    end

    for element in grid.indices[coord_indices...]
        if grid isa SusceptibilityGridNoRepeat || iszero(element.shift)
            shifted = normed
        else
            shifted = normed .- grid.shifts[element.shift]
        end
        field += element_susceptibility(element, grid, shifted, inside)
    end
    return field
end

function susceptibility_off_resonance(fixed::FixedSusceptibility, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing)
    sum(susceptibility_off_resonance(grid, position, inside) for grid in fixed)
end

susceptibility_off_resonance(fixed::FixedSusceptibility{0}, position::SVector{3, Float64}, inside::Union{Nothing, Bool}=nothing) = 0.


"""
    dipole_approximation(susceptibility, offset, distance, B0_field)
    dipole_approximation(magnetisation, offset, distance, B0_field)

Computes the shift in magnetic field due to a shift of susceptibility which is at the given offset (with given pre-computed distance) along the B0 field.

The `susceptibility` is given as a scalar value and is assumed to generate a field in the z-direction.
The `magnetisation` is a vector representing a magnetisation in an arbitrary direction (as required for anisotropic susceptibilities).

Equations from [schenck96_role_magnet_suscep_magnet_reson_imagin](@cite).
"""
function dipole_approximation(susceptiblity::Float64, offset::SVector{3, Float64}, dist::Float64, B0_field::SVector{3, Float64})
    if iszero(dist)
        return 0.
    end
    offset_B0 = (
        offset[1] * B0_field[1] +
        offset[2] * B0_field[2] +
        offset[3] * B0_field[3]
    )
    return susceptiblity / (4π * dist^5) * (3 * offset_B0^2 - dist^2)
end

function dipole_approximation(magnetisation::SVector{3, Float64}, offset::SVector{3, Float64}, dist::Float64, B0_field::SVector{3, Float64})
    offset_mag = (
        magnetisation[1] * offset[1] +
        magnetisation[2] * offset[2] +
        magnetisation[3] * offset[3]
    )
    offset_b0 = (
        B0_field[1] * offset[1] +
        B0_field[2] * offset[2] +
        B0_field[3] * offset[3]
    )
    mag_b0 = (
        B0_field[1] * magnetisation[1] +
        B0_field[2] * magnetisation[2] +
        B0_field[3] * magnetisation[3]
    )
    return 1 / (4π * dist^5) * (3 * offset_b0 * offset_mag - mag_b0 * dist^2)
end

average_field(x, y, nsteps) = mean([1/sqrt((d+x)^2 + y^2)^3 for d in range(-0.5, 0.5, nsteps * 2 + 1)[2:2:end-1]])

function dipole_approximation(susceptiblity::Float64, offset::SVector{2, Float64}, dist::Float64, B0_field::SVector{2, Float64})
    if iszero(dist)
        return 0.
    end
    offset_B0 = (
        offset[1] * B0_field[1] +
        offset[2] * B0_field[2]
    )
    inplane_B0_sqr = (
        B0_field[1] * B0_field[1] +
        B0_field[2] * B0_field[2]
    )
    if iszero(inplane_B0_sqr)
        return 0.
    end
    return susceptiblity / (2π * dist^4) * (2 * offset_B0^2 - dist^2 * inplane_B0_sqr)
end

function dipole_approximation(magnetisation::SVector{2, Float64}, offset::SVector{2, Float64}, dist::Float64, B0_field::SVector{2, Float64})
    offset_mag = (
        magnetisation[1] * offset[1] +
        magnetisation[2] * offset[2]
    )
    offset_b0 = (
        B0_field[1] * offset[1] +
        B0_field[2] * offset[2]
    )
    mag_b0 = (
        B0_field[1] * magnetisation[1] +
        B0_field[2] * magnetisation[2]
    )
    return 1 / (2π * dist^4) * (2 * offset_b0 * offset_mag - mag_b0 * dist^2)
end

function dipole_approximation_repeat(susceptiblity::Float64, offset::SVector{N, Float64}, B0_field::SVector{N, Float64}, repeats::SVector{N, Float64}) where{N}
    half_repeats = repeats ./ 2
    normed = @. mod(offset + half_repeats, repeats) - half_repeats
    result = 0.
    for shift in Iterators.product(Iterators.repeated(-1:1, N)...)
        shifted_offset = @. normed + shift * repeats
        result += dipole_approximation(susceptiblity, shifted_offset, norm(shifted_offset), B0_field)
    end
    return result
end

"""
    element_susceptibility(source::SuscetibilityGridElement, grid::SuscetibilityGrid, position, stuck_to)

Computes the off-resonance field contribution from a [`SuscetibilityGridElement`](@ref).

For a `position` within twice the `source.radius` of `source.position`, this will call [`dipole_approximation`](@ref).
For any closer `position` [`single_susceptibility`](@ref) will be called on the appropriate element in `grid.sources`.
"""
function element_susceptibility(element::SusceptibilityGridElement, grid::SusceptibilityGrid, position::AbstractVector, stuck_inside::Union{Nothing, Bool})
    offset = position - element.position
    dist = norm(offset)

    if dist > element.radius
        return dipole_approximation(element.susceptibility, offset, dist, grid.B0_field)
    else
        return single_susceptibility(
            grid.sources[element.index],
            offset,
            dist,
            stuck_inside,
            grid.B0_field
        )
    end
end

"""
    off_resonance_gradient(susceptibility, B0)

Maximum gradient of the off-resonance field in kHz/um due to the simulated susceptibility sources.

Internally, computed for each susceptibility sources using [`single_susceptibility_gradient`](@ref).
The maximum out of these is returned.
"""
function off_resonance_gradient(grid::SusceptibilityGrid, B0)
    return maximum(single_susceptibility_gradient.(grid.sources)) * B0 * gyromagnetic_ratio * 1e-6
end

function off_resonance_gradient(fixed::FixedSusceptibility, B0)
    maximum(off_resonance_gradient(grid, B0) for grid in fixed)
end

off_resonance_gradient(::FixedSusceptibility{0}, B0) = 0.



end