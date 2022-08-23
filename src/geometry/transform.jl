using StaticArrays
struct TransformObstruction{N, M, K, O<:BaseObstruction{N}} <: Obstruction{N}
    obstructions::SVector{M, O}
    shifts::SVector{M, SVector{N, Float}}
    repeats::SVector{N, Float}
    shift_quadrants::SVector{M, SVector{N, Bool}}
    rotation::SMatrix{3, N, Float, K}
    inv_rotation::SMatrix{N, 3, Float, K}
    chi::Float
    lorentz_radius :: Float
    lorentz_repeats :: SVector{N, Int}
    function TransformObstruction(obstructions::SVector{M, <:BaseObstruction{N}}, shifts, repeats, rotation, lorentz_radius) where {N, M}
        @assert all(r -> r>0, repeats)
        repeats = map(Float, repeats)

        # normalise shifts by the repeats, so that the base shift is between -repeats/4 and +3 * repeats/4
        shifts = SVector{M}([map((s, r) -> isfinite(r) ? mod(s + r/4, r) - r/4 : Float(s), single_shift, repeats) for single_shift in shifts])
        shift_quadrants = SVector{M}([get_shift_quadrants(s, repeats) for s in shifts])
        chi = sum(total_susceptibility, obstructions) / prod(repeats)
        lorentz_repeats = map(r -> isfinite(r) ? Int(div(lorentz_radius, r, RoundUp)) : 0, repeats)
        new{N, M, 3 * N, eltype(obstructions)}(obstructions, shifts, repeats, shift_quadrants, rotation, transpose(rotation), chi, lorentz_radius, lorentz_repeats)
    end
end

get_shift_quadrants(shift::SVector{N, Float}, repeats::SVector{N, Float}) where {N} = map((s, r) -> isfinite(r) && s > r/4., shift, repeats)

function TransformObstruction(Obstruction::Type{<:BaseObstruction{N}}, args...; shifts=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10., kwargs...) where {N}
    o = Obstruction.(args...; kwargs...)
    TransformObstruction(o; shifts=shifts, repeats=repeats, rotation=rotation, lorentz_radius=lorentz_radius)
end

function TransformObstruction(obstruction::BaseObstruction{N}; shifts=nothing, kwargs...) where N
    if isnothing(shifts)
        return TransformObstruction([obstruction]; shifts=shifts, kwargs...)
    elseif isa(shifts, Real)
        @assert N == 1
        return TransformObstruction([obstruction]; shifts=[shifts], kwargs...)
    elseif eltype(shifts) <: Real && N > 1
        @assert length(shifts) == N
        return TransformObstruction([obstruction]; shifts=[shifts], kwargs...)
    else
        return TransformObstruction([copy(obstruction) for _ in 1:length(shifts)]; shifts=shifts, kwargs...)
    end
end

get_rotation(rotation::Rotations.Rotation, ndim::Int) = get_rotation(Rotations.RotMatrix(rotation.mat), ndim)
get_rotation(rotation::Rotations.RotMatrix, ndim::Int) = get_rotation(rotation.mat, ndim)
function get_rotation(rotation::AbstractMatrix, ndim::Int)
    if size(rotation) == (3, 3) && ndim < 3
        rotation = rotation[:, 1:ndim]
    end
    SMatrix{3, ndim, Float}(rotation)
end

function get_rotation(rotation::AbstractVector, ndim::Int)
    @assert length(rotation)==3
    normed = rotation / norm(rotation)
    if ndim == 1
        return reshape(normed, 3, 1)
    end
    try_vec = [1., 0, 0]
    vec1 = cross(normed, try_vec)
    if iszero(norm(vec1))
        try_vec = [0., 1, 0]
        vec1 = cross(normed, try_vec)
    end
    vec2 = cross(normed, vec1)
    vec1 = vec1 ./ norm(vec1)
    vec2 = vec2 ./ norm(vec2)
    if ndim == 2
        return hcat(vec1, vec2)
    else
        return hcat(vec1, vec2, normed)
    end
end

function get_rotation(rotation::Symbol, ndim::Int)
    target_dimension = Dict(
        :x => 1,
        :y => 2,
        :z => 3,
    )[rotation]
    orig_dimension = ndim == 1 ? 1 : 3
    target = zeros(Float, 3, ndim)
    if orig_dimension <= ndim
        target[target_dimension, orig_dimension] = one(Float)
    end
    if target_dimension <= ndim
        target[orig_dimension, target_dimension] = one(Float)
    end
    for d in 1:ndim
        if d != target_dimension && d != orig_dimension
            target[d, d] = one(Float)
        end
    end
    SMatrix{3, ndim, Float}(target)
end

function TransformObstruction(obstructions::AbstractVector{<:BaseObstruction{N}}; shifts=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10) where {N}
    obstructions = Vector(obstructions)
    if isnothing(shifts)
        shifts = zeros(SVector{N, Float}, length(obstructions))
    else
        if eltype(shifts) <: Real && N > 1
            error("Shifts cannot by scalar numbers for $(eltype(obstructions))")
        end
        shifts = [SVector{N, Float}(a) for a in shifts]
    end
    if isa(repeats, Real)
        repeats = fill(repeats, N)
    end
    M = length(obstructions)
    TransformObstruction(
        SVector{M}(obstructions),
        SVector{M}(shifts),
        SVector{N}([iszero(r) ? Float(Inf) : Float(r) for r in repeats]),
        get_rotation(rotation, N),
        lorentz_radius
    )
end

function project_rotation(trans::TransformObstruction{N}, pos::PosVector) where {N}
    trans.inv_rotation * pos
end

function project_repeat(trans::TransformObstruction{N}, pos::SVector{N, Float}) where {N}
    pos_shifted = map((p, r) -> isfinite(r) ? mod(p, r) : p, pos, trans.repeats)
    pos_center = map((p, r) -> (isfinite(r) && p > r/2) ? p - r : p, pos_shifted, trans.repeats)
    [map(project_repeat_scalar, pos_shifted, pos_center, s, q) for (s, q) in zip(trans.shifts, trans.shift_quadrants)]
end

function project_repeat_scalar(pos_shifted::Float, pos_center:: Float, shift::Float, shift_quadrant::Bool) where {N}
    (shift_quadrant ? pos_shifted : pos_center) - shift
end

function project(trans::TransformObstruction{N}, pos::PosVector) where {N}
    project_repeat(trans, project_rotation(trans, pos))
end

function detect_collision(movement :: Movement{3}, trans :: TransformObstruction{N}, previous=empty_collision) where {N}
    # apply rotation
    projected_origin = project_rotation(trans, movement.origin)
    projected_destination = project_rotation(trans, movement.destination)
    if !any(isfinite, trans.repeats)
        c = empty_collision

        for (shift, obstruction) in zip(trans.shifts, trans.obstructions)
            ctest = detect_collision(
                Movement{N}(projected_origin - shift, projected_destination - shift, Float(1)),
                obstruction,
                previous
            )
            if ctest.distance < c.distance
                c = ctest
            end
        end
    elseif all(isfinite, trans.repeats)
        c = detect_collision(projected_origin, projected_destination, previous, trans.repeats, trans.shifts, trans.shift_quadrants, trans.obstructions)
    else
        error()
    end
    if c.distance <= 1.
        if N == 1
            n = SVector{1, Float}(c.normal[1])
        elseif N == 2
            n = SVector{2, Float}(c.normal[1], c.normal[2])
        else
            n = c.normal
        end
        return Collision(
            c.distance,
            trans.rotation * n,
            c.id,
            c.index
        )
    end
    return empty_collision
end

    
function detect_collision(
    origin :: SVector{N, Float}, destination :: SVector{N, Float}, previous::Collision, repeats::SVector{N, Float}, 
    shifts::SVector{M, SVector{N, Float}}, shift_quadrants::SVector{M, SVector{N, Bool}},
    obstructions::SVector{M, <:BaseObstruction}
    )  where {N, M}
    half_repeats = map(r -> 0.5 * r, repeats)
    # project onto 1x1x1 grid
    grid_origin = map(/, origin, half_repeats)
    grid_destination = map(/, destination, half_repeats)

    # iterate over multiple crossings of repeat boundaries
    for (voxel, _, _, t2, _) in ray_grid_intersections(grid_origin, grid_destination)
        lower = map(v -> mod(v, 2), voxel)
        
        voxel_orig = map(*, voxel + lower, half_repeats)
        p1 = origin - voxel_orig
        p2 = destination - voxel_orig

        c = empty_collision
        for (shift, quadrant, obstruction) in zip(shifts, shift_quadrants, obstructions)
            to_shift = map((s, q, h) -> s + q * h, shift, quadrant, half_repeats)
            ctest = detect_collision(
                Movement{N}(p1 - to_shift, p2 - to_shift, Float(1)),
                obstruction,
                previous
            )
            if ctest.distance < c.distance
                c = ctest
            end
        end
        if c.distance <= t2
            return c
        end
    end
    return empty_collision
end


isinside(trans::TransformObstruction, pos::PosVector) = any(po -> isinside(po[1], po[2]), zip(project(trans, pos), trans.obstructions))

function BoundingBox(trans::TransformObstruction{N}) where {N}
    bbs = [BoundingBox(o) for o in trans.obstructions]
    lower_bounds = [bb.lower + shift for (bb, shift) in zip(bbs, trans.shifts)]
    upper_bounds = [bb.upper + shift for (bb, shift) in zip(bbs, trans.shifts)]
    bb = BoundingBox(
        SVector{N, Float}([minimum(b -> b[dim], lower_bounds) for dim in 1:N]),
        SVector{N, Float}([maximum(b -> b[dim], upper_bounds) for dim in 1:N]),
    )
end


function off_resonance(transform::TransformObstruction{N}, position::PosVector, b0_field::PosVector=PosVector([0, 0, 1])) where {N}
    b0 = project_rotation(transform, b0_field)

    # Contribution from outside of the Lorentz cavity
    if N == 2
        sinsq = b0[1] * b0[1] + b0[2] * b0[2]
        field = transform.chi * sinsq
    else
        field = transform.chi
    end

    finite_repeats = map(r -> isfinite(r) ? r : zero(Float), transform.repeats)
    # Contribution from within the Lorentz cavity
    return field + sum(map((o, p)->lorentz_off_resonance(o, p, b0, finite_repeats, transform.lorentz_radius, transform.lorentz_repeats), transform.obstructions, project(transform, position)))
end