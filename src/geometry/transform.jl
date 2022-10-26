using StaticArrays
"""
    TransformObstruction(base_type, args...; positions=nothing, repeats=0., rotation=I(3), lorentz_radius=10., kwargs...)
    TransformObstruction(base; positions=nothing, repeats=0., rotation=I(3), lorentz_radius=10.)

Transforms the [`BaseObstruction`](@ref) objects in `base` in one of several ways:
- Apply the shifts defined by `positions` (none applied by default).
- Infinitely repeat the base obstruction every `repeats`.
- Rotate the base obstructions by `rotation`.

The [`BaseObstruction`](@ref) objects an already be initialised.
If a [`BaseObstruction`](@ref) type is given instead it will be initialised using `args` and `kwargs` not used by `TransformObstruction`.
"""
struct TransformObstruction{N, M, K, O<:BaseObstruction{N}} <: Obstruction{N}
    obstructions::NTuple{M, O}
    positions::NTuple{M, SVector{N, Float}}
    repeats::SVector{N, Float}
    shift_quadrants::NTuple{M, SVector{N, Bool}}
    rotation::SMatrix{3, N, Float, K}
    inv_rotation::SMatrix{N, 3, Float, K}
    chi::Float
    lorentz_radius :: Float
    lorentz_repeats :: SVector{N, Int}
    function TransformObstruction(obstructions::NTuple{M, <:BaseObstruction{N}}, positions, repeats, rotation, lorentz_radius) where {N, M}
        @assert all(r -> r>0, repeats)
        repeats = SVector{N}(map(Float, repeats))

        # normalise shifts by the repeats, so that the base shift is between -repeats/4 and +3 * repeats/4
        positions = NTuple{M, SVector{N, Float}}([map((s, r) -> isfinite(r) ? mod(s + r/4, r) - r/4 : Float(s), single_shift, repeats) for single_shift in positions])
        shift_quadrants = NTuple{M, SVector{N, Bool}}([get_shift_quadrants(s, repeats) for s in positions])
        chi = sum(total_susceptibility, obstructions) / prod(repeats)
        lorentz_repeats = map(r -> isfinite(r) ? Int(div(lorentz_radius, r, RoundUp)) : 0, repeats)
        new{N, M, 3 * N, eltype(obstructions)}(obstructions, positions, repeats, shift_quadrants, rotation, transpose(rotation), chi, lorentz_radius, lorentz_repeats)
    end
end

get_shift_quadrants(shift::SVector{N, Float}, repeats::SVector{N, Float}) where {N} = map((s, r) -> isfinite(r) && s > r/4., shift, repeats)

function TransformObstruction(Obstruction::Type{<:BaseObstruction{N}}, args...; positions=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10., kwargs...) where {N}
    o = Obstruction.(args...; kwargs...)
    TransformObstruction(o; positions=positions, repeats=repeats, rotation=rotation, lorentz_radius=lorentz_radius)
end

function TransformObstruction(obstruction::BaseObstruction{N}; positions=nothing, kwargs...) where N
    if isnothing(positions)
        return TransformObstruction([obstruction]; positions=positions, kwargs...)
    elseif isa(positions, Real)
        @assert N == 1
        return TransformObstruction([obstruction]; positions=[positions], kwargs...)
    elseif eltype(positions) <: Real && N > 1
        @assert length(positions) == N
        return TransformObstruction([obstruction]; positions=[positions], kwargs...)
    else
        return TransformObstruction([copy(obstruction) for _ in 1:length(positions)]; positions=positions, kwargs...)
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
    try_vec = [0., 1., 0.]
    vec1 = cross(normed, try_vec)
    if iszero(norm(vec1))
        try_vec = [1., 0, 0]
        vec1 = cross(normed, try_vec)
    end
    vec2 = cross(normed, vec1)
    vec1 = vec1 ./ norm(vec1)
    vec2 = vec2 ./ norm(vec2)
    return get_rotation(hcat(vec1, vec2, normed), ndim)
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

function TransformObstruction(obstructions::AbstractVector{<:BaseObstruction{N}}; positions=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10) where {N}
    obstructions = Vector(obstructions)
    if isnothing(positions)
        positions = zeros(SVector{N, Float}, length(obstructions))
    else
        if eltype(positions) <: Real && N > 1
            error("Positions cannot by scalar numbers for $(eltype(obstructions))")
        end
        positions = [SVector{N, Float}(a) for a in positions]
    end
    if isa(repeats, Real)
        repeats = fill(repeats, N)
    end
    M = length(obstructions)
    TransformObstruction(
        NTuple{M, eltype(obstructions)}(obstructions),
        NTuple{M, SVector{N, Float}}(positions),
        SVector{N, Float}([iszero(r) ? Float(Inf) : Float(r) for r in repeats]),
        get_rotation(rotation, N),
        lorentz_radius
    )
end

function project_rotation(trans::TransformObstruction{N}, pos::PosVector) where {N}
    trans.inv_rotation * pos
end

function project_repeat(trans::TransformObstruction{N, M}, pos::SVector{N, Float}) where {N, M}
    pos_shifted = map((p, r) -> isfinite(r) ? mod(p, r) : p, pos, trans.repeats)
    will_shift = map((p, r) -> (isfinite(r) && p > r/2) ? r : zero(Float), pos_shifted, trans.repeats)

    function shifted(quadrant, shift)
        toshift = map(+, map((q, r) -> (!q) * r, quadrant, will_shift), shift)
        return map(-, pos_shifted, toshift)
    end
    map(shifted, trans.shift_quadrants, trans.positions) :: NTuple{M, SVector{2, Float}}
end

function project(trans::TransformObstruction{N}, pos::PosVector) where {N}
    project_repeat(trans, project_rotation(trans, pos))
end

function detect_collision(movement :: Movement{3}, trans :: TransformObstruction{N}, previous::Collision) where {N}
    # apply rotation
    projected_origin = project_rotation(trans, movement.origin)
    projected_destination = project_rotation(trans, movement.destination)
    if !any(isfinite, trans.repeats)
        c = empty_collision

        for (shift, obstruction) in zip(trans.positions, trans.obstructions)
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
        c = detect_collision(projected_origin, projected_destination, previous, trans.repeats, trans.positions, trans.shift_quadrants, trans.obstructions)
    else
        error("Repeating only along a subset of directions is not yet supported")
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
            c.properties,
            c.index
        )
    end
    return empty_collision
end

    
function detect_collision(
    origin :: SVector{N, Float}, destination :: SVector{N, Float},
    previous::Collision, repeats::SVector{N, Float}, 
    positions::NTuple{M, SVector{N, Float}}, shift_quadrants::NTuple{M, SVector{N, Bool}},
    obstructions::NTuple{M, <:BaseObstruction}
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
        to_correct = map(*, lower, repeats)

        c = empty_collision
        for (shift, quadrant, obstruction) in zip(positions, shift_quadrants, obstructions)
            to_shift = map((s, q, t) -> s - q * t, shift, quadrant, to_correct)
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


isinside(trans::TransformObstruction, pos::PosVector) = maximum(po -> isinside(po[2], po[1]), zip(project(trans, pos), trans.obstructions))

function BoundingBox(trans::TransformObstruction{N}) where {N}
    bbs = [BoundingBox(o) for o in trans.obstructions]
    lower_bounds = [bb.lower + shift for (bb, shift) in zip(bbs, trans.positions)]
    upper_bounds = [bb.upper + shift for (bb, shift) in zip(bbs, trans.positions)]
    bb = BoundingBox(
        SVector{N, Float}([minimum(b -> b[dim], lower_bounds) for dim in 1:N]),
        SVector{N, Float}([maximum(b -> b[dim], upper_bounds) for dim in 1:N]),
    )
end


function off_resonance(transform::TransformObstruction{N, M}, position::PosVector, b0_field::PosVector=PosVector([0, 0, 1])) where {N, M}
    b0 = project_rotation(transform, b0_field)

    finite_repeats = map(r -> isfinite(r) ? r : zero(Float), transform.repeats)
    # Contribution from within the Lorentz cavity
    off_resonances = map((o, p)->lorentz_off_resonance(o, p, b0, finite_repeats, transform.lorentz_radius, transform.lorentz_repeats), transform.obstructions, project(transform, position)) :: NTuple{M, Float}
    return sum(off_resonances)
end