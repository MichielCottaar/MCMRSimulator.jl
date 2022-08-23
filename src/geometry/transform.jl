struct TransformObstruction{N, O<:BaseObstruction{N}} <: Obstruction{N}
    obstructions::Vector{O}
    shifts::Vector{SVector{N, Float}}
    repeats::SVector{N, Float}
    shift_quadrants::Vector{SVector{N, Bool}}
    rotation::SMatrix{3, N, Float}
    inv_rotation::SMatrix{N, 3, Float}
    chi::Float
    lorentz_radius :: Float
    lorentz_repeats :: SVector{N, Int}
    function TransformObstruction(obstructions::Vector{<:BaseObstruction{N}}, shifts, repeats, rotation, lorentz_radius) where {N}
        @assert all(r -> r>0, repeats)
        repeats = map(Float, repeats)

        # normalise shifts by the repeats, so that the base shift is between -repeats/4 and +3 * repeats/4
        shifts = [map((s, r) -> isfinite(r) ? mod(s + r/4, r) - r/4 : Float(s), single_shift, repeats) for single_shift in shifts]
        shift_quadrants = [get_shift_quadrants(s, repeats) for s in shifts]
        chi = sum(total_susceptibility, obstructions) / prod(repeats)
        lorentz_repeats = map(r -> isfinite(r) ? Int(div(lorentz_radius, r, RoundUp)) : 0, repeats)
        new{N, eltype(obstructions)}(obstructions, shifts, repeats, shift_quadrants, rotation, transpose(rotation), chi, lorentz_radius, lorentz_repeats)
    end
end

get_shift_quadrants(shift::SVector{N, Float}, repeats::SVector{N, Float}) where {N} = map((s, r) -> isfinite(r) && s > r/4., shift, repeats)

function TransformObstruction(Obstruction::Type{<:BaseObstruction{N}}, args...; shifts=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10., kwargs...) where {N}
    o = Obstruction.(args...; kwargs...)
    TransformObstruction(o; shifts=shifts, repeats=repeats, rotation=rotation, lorentz_radius=lorentz_radius)
end

function TransformObstruction(obstruction::BaseObstruction{N}; shifts, kwargs...) where N
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
    TransformObstruction(
        obstructions,
        shifts,
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
    
    # project onto 1x1x1 grid
    fdiv(p, r) = isfinite(r) ? 2 * p / r : 0.5
    frev(r, p, lower, p_orig) = isfinite(r) ? r * (p - lower) / 2 : p_orig
    origin = map(fdiv, projected_origin, trans.repeats)
    destination = map(fdiv, projected_destination, trans.repeats)

    # iterate over multiple crossings of repeat boundaries
    for (voxel, t1, p1, t2, p2) in ray_grid_intersections(origin, destination)
        lower = mod.(voxel, 2)
        upper = 1 .- lower
        orig1 = (projected_destination .* t1) .+ (projected_origin .* (1 - t1))
        orig2 = (projected_destination .* t2) .+ (projected_origin .* (1 - t2))
        pos1_center = map(frev, trans.repeats, p1, lower, orig1)
        pos2_center = map(frev, trans.repeats, p2, lower, orig2)
        pos1_shift = map(frev, trans.repeats, p1, upper, orig1)
        pos2_shift = map(frev, trans.repeats, p2, upper, orig2)

        c = empty_collision
        for (s, q, o) in zip(trans.shifts, trans.shift_quadrants, trans.obstructions)
            part_origin = map(project_repeat_scalar, pos1_shift, pos1_center, s, q)
            part_destination = map(project_repeat_scalar, pos2_shift, pos2_center, s, q)
            ctest = detect_collision(
                Movement{N}(part_origin, part_destination, Float(1)),
                o,
                previous
            )
            if ctest.distance < c.distance
                c = ctest
            end
        end
        if c !== empty_collision
            return Collision(
                c.distance * (t2 - t1) + t1,
                trans.rotation * SVector{N, Float}(c.normal[1:N]),
                c.id,
                c.index
            )
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