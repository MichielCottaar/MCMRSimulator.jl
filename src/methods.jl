"""
Defines methods shared across multiple sub-modules.
"""
module Methods

import Rotations
import LinearAlgebra: norm, cross, ⋅
import StaticArrays: SMatrix

"""
    get_time(snapshot)
    get_time(sequence_component)
    get_time(sequence, sequence_index)

Returns the time in milliseconds that a snapshot was taken or that a sequence component will have effect.
"""
function get_time end

"""
    B0(scanner)
    B0(sequence)

Returns the magnetic field strength of the scanner in Tesla.
"""
function B0 end

# Used for spins and RF pulses
function phase end

# Used for plot_plane, spin, and geometry
function project end

"""
    norm_angle(angle)

Normalises an angle in degrees, so that it is between it is in the range (-180, 180]
"""
function norm_angle(angle)
    angle = mod(angle, 360)
    if angle > 180
        angle -= 360
    end
    angle
end


function off_resonance end

"""
    get_rotation(rotation_mat, ndim)

Returns a (3, `ndim`) rotation matrix, that is the relevant part of the full 3x3 `rotation_mat` to map to the x-axis (if `ndim` is 1) or the x-y plance (if `ndim` is 2).
If `ndim` is 3, the full rotation matrix `rotation_mat` is returned.

    get_rotation(vector, ndim)

Returns the (3, `ndim`) rotation matrix mapping the `vector` to the x-direction (if `ndim` is 1 or 3) or the z-direction (if `ndim` is 2).
Vector can be a length 3 array or one of the symbols :x, :y, or :z (representing vectors in those cardinal directions).
"""
get_rotation(rotation::Rotations.Rotation, ndim::Int) = get_rotation(Rotations.RotMatrix(rotation), ndim)
get_rotation(rotation::Rotations.RotMatrix, ndim::Int) = get_rotation(rotation.mat, ndim)
function get_rotation(rotation::AbstractMatrix, ndim::Int)
    if size(rotation) == (3, 3) && ndim < 3
        rotation = rotation[:, 1:ndim]
    end
    SMatrix{3, ndim, Float64}(rotation)
end

function get_rotation(rotation::AbstractVector, ndim::Int)
    T = typeof(rotation[1])
    return get_rotation(T.(rotation), ndim)
end

function get_rotation(rotation::AbstractVector{<:AbstractVector}, ndim::Int)
    return get_rotation(hcat(rotation...), ndim)
end

function get_rotation(rotation::AbstractVector{<:Number}, ndim::Int)
    @assert length(rotation)==3
    normed = rotation / norm(rotation)
    if ndim == 1
        return reshape(normed, 3, 1)
    end
    try_vec = [0., 1., 0.]
    vec1 = cross(normed, try_vec)
    if iszero(norm(vec1))
        try_vec = [1., 0, 0]
        vec1 = cross(normed, try_vec)
    end
    reference = ndim == 2 ? [0, 0, 1.] : [1., 0., 0.]
    # rotate reference to normed
    angle = acos(normed ⋅ reference)
    vec_rotation = cross(reference, normed)
    if iszero(norm(vec_rotation))
        return get_rotation(abs(angle) < 1 ? I(3) : -I(3), 3)
    end
    return get_rotation(Rotations.AngleAxis(angle, vec_rotation...), ndim)
end

function get_rotation(rotation::Symbol, ndim::Int)
    orig_dimension = ndim == 1 ? 1 : 3
    target_dimension = Dict(
        :x => 1,
        :y => 2,
        :z => 3,
        :I => orig_dimension,
    )[rotation]
    target = zeros(Float64, 3, ndim)
    if orig_dimension <= ndim
        target[target_dimension, orig_dimension] = one(Float64)
    end
    if target_dimension <= ndim
        target[orig_dimension, target_dimension] = one(Float64)
    end
    for d in 1:ndim
        if d != target_dimension && d != orig_dimension
            target[d, d] = one(Float64)
        end
    end
    SMatrix{3, ndim, Float64}(target)
end

end