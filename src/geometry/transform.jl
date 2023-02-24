using StaticArrays
struct SingleTransform{N, B<:BoundingBox{N}}
    shift::SVector{N, Float}
    quadrant::SVector{N, Bool}
    bounding_box::B
end

"""
    TransformObstruction(base; positions=nothing, repeats=0., rotation=I(3), lorentz_radius=10.)
    TransformObstruction(base_type, args...; positions=nothing, repeats=0., rotation=I(3), lorentz_radius=10., kwargs...)

Transforms the [`BaseObstruction`](@ref) objects in `base` in one of several ways:
- Apply the shifts defined by `positions` (none applied by default).
- Infinitely repeat the base obstruction every `repeats`.
- Rotate the base obstructions by `rotation`.

If a [`BaseObstruction`](@ref) type is given instead it will be initialised using any `args` and `kwargs` not used by `TransformObstruction`.
"""
struct TransformObstruction{N, M, K, B<:BoundingBox{N}, O<:BaseObstruction{N}} <: Obstruction{N}
    transforms::Vector{SingleTransform{N, B}}
    obstructions::Vector{O}
    repeats::SVector{N, Float}
    rotation::SMatrix{3, N, Float, K}
    inv_rotation::SMatrix{N, 3, Float, K}
    chi::Float
    lorentz_radius :: Float
    lorentz_repeats :: SVector{N, Int}
    function TransformObstruction(obstructions::AbstractVector{<:BaseObstruction{N}}, positions, repeats, rotation, lorentz_radius) where {N}
        @assert all(r -> r>0, repeats)
        repeats = SVector{N}(map(Float, repeats))

        # normalise shifts by the repeats, so that the base shift is between -repeats/4 and +3 * repeats/4
        positions = Vector([map((s, r) -> isfinite(r) ? mod(s + r/4, r) - r/4 : Float(s), single_shift, repeats) for single_shift in positions])
        shift_quadrants = [get_shift_quadrants(s, repeats) for s in positions]
        chi = sum(total_susceptibility, obstructions) / prod(repeats)
        lorentz_repeats = map(r -> isfinite(r) ? Int(div(lorentz_radius, r, RoundUp)) : 0, repeats)
        @assert length(obstructions) == length(positions)
        bounding_boxes = BoundingBox.(obstructions)
        for (pos, bb, qq) in zip(positions, bounding_boxes, shift_quadrants)
            u = upper(bb)
            l = lower(bb)
            for dim in 1:N
                upper_limit = qq[dim] ? repeats[dim] : (repeats[dim] / 2)
                lower_limit = qq[dim] ? zero(Float) : -(repeats[dim] / 2)
                if (
                    (u[dim] + pos[dim]) > upper_limit ||
                    (l[dim] + pos[dim]) < lower_limit
                )
                    throw(DomainError("$(eltype(obstructions)) object at position $(pos) is too large along dimension $(dim) to fit within repeats of $(repeats[dim])um."))
                end
            end
        end
        transforms = SingleTransform{N, eltype(bounding_boxes)}.(positions, shift_quadrants, bounding_boxes)
        new{N, length(obstructions), 3 * N, eltype(bounding_boxes), eltype(obstructions)}(transforms, obstructions, repeats, rotation, transpose(rotation), chi, lorentz_radius, lorentz_repeats)
    end
end

function Base.show(io::IO, geom::TransformObstruction{N, M, K, B, O}) where {N, M, K, B, O}
    print(io, "$(length(geom.obstructions)) ")
    if any(isfinite.(geom.repeats))
        print(io, "repeating ")
    end
    print(io, String(nameof(O)) * " objects")
end
get_shift_quadrants(shift::SVector{N, Float}, repeats::SVector{N, Float}) where {N} = map((s, r) -> isfinite(r) && s > r/4., shift, repeats)

function TransformObstruction(Obstruction::Type{<:BaseObstruction{N}}, args...; positions=nothing, repeats=0., rotation=one(Rotations.RotMatrix{3}), lorentz_radius=10., kwargs...) where {N}
    n_args = length(args)
    keywords = keys(kwargs)
    function call_constructor(arguments...)
        args = arguments[1:n_args]
        keyword_values = arguments[n_args+1:end]
        kwargs = Dict(symbol => value for (symbol, value) in zip(keywords, keyword_values))
        Obstruction(args...; kwargs...)
    end
    o = call_constructor.(args..., values(kwargs)...)
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
        return TransformObstruction([deepcopy(obstruction) for _ in 1:length(positions)]; positions=positions, kwargs...)
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

full_rotation_mat(transform::TransformObstruction{3}) = transform.rotation
function full_rotation_mat(transform::TransformObstruction{2})
    v1 = transform.rotation[:, 1]
    v2 = transform.rotation[:, 2]
    v3 = cross(v1, v2)
    return SMatrix{3, 3, Float, 9}(hcat(v1, v2, v3))
end
function full_rotation_mat(transform::TransformObstruction{1})
    v1 = transform.rotation[:, 1]
    v2 = v1
    for try_v2 in ([1, 0, 0], [0, 1, 0])
        v2_unnorm = cross(v1, try_v2)
        v2_norm = norm(v2_unnorm)
        if iszero(v2_norm)
            continue
        end
        v2 = v2_unnorm / v2_norm
    end
    v3 = cross(v1, v2)
    return SMatrix{3, 3, Float, 9}(hcat(v1, v2, v3))
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
    TransformObstruction(
        obstructions,
        positions,
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

    function shifted(single::SingleTransform)
        toshift = map(+, map((q, r) -> (!q) * r, single.quadrant, will_shift), single.shift)
        return map(-, pos_shifted, toshift)
    end
    return shifted
end

function project(trans::TransformObstruction{N, M}, pos::PosVector) where {N, M}
    [project_repeat(trans, project_rotation(trans, pos))(trans.obstructions[index]) for index in 1:M]
end

function detect_collision(movement :: Movement{3}, trans :: TransformObstruction{N}, previous::Collision) where {N}
    # apply rotation
    projected_origin = project_rotation(trans, movement.origin)
    projected_destination = project_rotation(trans, movement.destination)
    if !any(isfinite, trans.repeats)
        current_guess = empty_collision

        for (index, single) in enumerate(trans.transforms)
            shifted_origin = projected_origin - single.shift
            shifted_destination = projected_destination - single.shift
            if possible_intersection(single.bounding_box, shifted_origin, shifted_destination)
                ctest = detect_collision(
                    Movement{N}(shifted_origin, shifted_destination, Float(1)),
                    trans.obstructions[index],
                    previous
                )
                if ctest.distance < current_guess.distance
                    current_guess = ctest
                end
            end
        end
    elseif all(isfinite, trans.repeats)
        current_guess = detect_collision_repeats(projected_origin, projected_destination, trans, previous)
    else
        error("Repeating only along a subset of directions is not yet supported")
    end
    if current_guess.distance <= 1.
        if N == 1
            n = SVector{1, Float}(current_guess.normal[1])
        elseif N == 2
            n = SVector{2, Float}(current_guess.normal[1], current_guess.normal[2])
        else
            n = current_guess.normal
        end
        return Collision(
            current_guess.distance,
            trans.rotation * n,
            current_guess.properties,
            current_guess.index,
            current_guess.inside
        )
    end
    return empty_collision
end

    
function detect_collision_repeats(
    origin :: SVector{N, Float}, destination :: SVector{N, Float},
    trans :: TransformObstruction{N, M},
    previous::Collision
    )  where {N, M}
    half_repeats = map(r -> 0.5 * r, trans.repeats)
    # project onto 1x1x1 grid
    grid_origin = map(/, origin, half_repeats)
    grid_destination = map(/, destination, half_repeats)
    if all(floor.(grid_origin) .== floor.(grid_destination))
        voxel = map(i->Int(floor(i)), grid_origin)
        return _detect_collision_repeats_helper(trans, voxel, half_repeats, origin, destination, previous)
    end

    # iterate over multiple crossings of repeat boundaries
    for (voxel, _, _, t2, _) in ray_grid_intersections(grid_origin, grid_destination)
        current_guess = _detect_collision_repeats_helper(trans, voxel, half_repeats, origin, destination, previous)
        if current_guess.distance <= t2
            return current_guess
        end
    end
    return empty_collision
end

@inline function _detect_collision_repeats_helper(trans::TransformObstruction{N, M}, voxel::SVector{N, Int}, half_repeats::SVector{N, Float}, origin::SVector{N, Float}, destination::SVector{N, Float}, previous::Collision) where {N, M}
    lower = map(v -> mod(v, 2), voxel)
        
    voxel_orig = map(*, voxel + lower, half_repeats)
    p1 = origin - voxel_orig
    p2 = destination - voxel_orig
    to_correct = map(*, lower, trans.repeats)

    current_guess = empty_collision
    for index in 1:M
        t = trans.transforms[index]
        to_shift = map((s, q, t) -> s - q * t, t.shift, t.quadrant, to_correct)
        shift1 = p1 - to_shift
        shift2 = p2 - to_shift
        if possible_intersection(t.bounding_box, shift1, shift2)
            ctest = detect_collision(
                Movement{N}(p1 - to_shift, p2 - to_shift, Float(1)),
                trans.obstructions[index],
                previous
            )
            if ctest.distance < current_guess.distance
                current_guess = ctest
            end
        end
    end
    return current_guess
end

function isinside(trans::TransformObstruction, pos::PosVector, stuck_to::Collision=empty_collision) 
    pos_func = project_repeat(trans, project_rotation(trans, pos))
    inside = 0
    for (index, t) in enumerate(trans.transforms)
        p = pos_func(t)
        if isinside(t.bounding_box, p)
            inside += isinside(trans.obstructions[index], p, stuck_to)
        end
    end
    return inside
end

function BoundingBox(trans::TransformObstruction{N}) where {N}
    lower_bounds = [lower(t.bounding_box) + t.shift for t in trans.transforms]
    upper_bounds = [upper(t.bounding_box) + t.shift for t in trans.transforms]
    return BoundingBox(
        SVector{N, Float}([minimum(b -> b[dim], lower_bounds) for dim in 1:N]),
        SVector{N, Float}([maximum(b -> b[dim], upper_bounds) for dim in 1:N]),
    )
end


function off_resonance(transform::TransformObstruction{N, M}, position::PosVector, b0_field::PosVector=PosVector([0, 0, 1])) where {N, M}
    b0 = project_rotation(transform, b0_field)

    finite_repeats = map(r -> isfinite(r) ? r : zero(Float), transform.repeats)

    # Contribution from within the Lorentz cavity
    positions_func = project_repeat(transform, project_rotation(transform, position))
    total = zero(Float)
    for index in 1:M
        total += lorentz_off_resonance(
            transform.obstructions[index],
            positions_func(transform.transforms[index]),
            b0, finite_repeats, transform.lorentz_radius,
            transform.lorentz_repeats
        )
    end
    return total
end

function produces_off_resonance(transform::TransformObstruction)
    for o in transform.obstructions
        if produces_off_resonance(o)
            return true
        end
    end
    return false
end


function inside_MRI_properties(transform::TransformObstruction, position::PosVector)
    positions_func = project_repeat(transform, project_rotation(transform, position))
    return merge_mri_parameters((
        inside_MRI_properties(obstruction, positions_func(single)) 
        for (obstruction, single) in zip(transform.obstructions, transform.transforms)
    ))
end


size_scale(t::TransformObstruction) = minimum(size_scale.(t.obstructions))
function size_scale(t::TransformObstruction{N, M, K, B, Wall}) where {N, M, K, B<:MCMRSimulator.BoundingBox{N}}
    min_dist = t.repeats[1]
    for s1 in t.transforms
        pos1 = s1.shift[1]
        for s2 in t.transforms
            pos2 = s2.shift[1]
            dist = abs(pos1 - pos2)
            if iszero(dist)
                continue
            end
            if dist > (t.repeats[1] / 2)
                dist = t.repeats[1] - dist
            end
            if dist < min_dist
                min_dist = dist
            end
        end
    end
    return min_dist
end

off_resonance_gradient(t::TransformObstruction) = maximum(off_resonance_gradient.(t.obstructions))
max_timestep_sticking(container::TransformObstruction, default_properties::GlobalProperties, diffusivity::Number) = maximum([max_timestep_sticking(o, default_properties, diffusivity) for o in container.obstructions])

function random_surface_spins(transform::TransformObstruction{N}, bounding_box::BoundingBox, volume_density::Number, default_surface_density; nsequences=1, kwargs...) where {N}
    rot_mat = full_rotation_mat(transform)
    inv_rot_mat = inv(rot_mat)
    rotated_bb = rotate(bounding_box, inv_rot_mat)
    all_spins = vcat([random_surface_spins(o, rotated_bb, volume_density, transform.repeats, default_surface_density, s.shift; kwargs...) for (o, s) in zip(transform.obstructions, transform.transforms)]...)
    for spin in all_spins
        spin.position = rot_mat * spin.position
    end
    return filter(spin->isinside(bounding_box, spin), all_spins)
end