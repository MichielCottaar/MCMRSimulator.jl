"""
Pseudo-volume sampling for mesh geometries.

Provides `sample_points_inside` for generating point distributions on mesh surfaces
and within mesh volumes using blue noise (Poisson disk) sampling.
"""
module Sampling

import Random
import LinearAlgebra: norm, ×
import StaticArrays: SVector
import Distributions: Poisson
import NearestNeighbors: KDTree, knn
import ...Geometries
import ...Internal: FixedGeometry, isinside, size_scale, geometry_mesh, PhysicalGeometry
import ...Internal.InternalBoundingBoxes: lower, upper
import ...Internal.PhysicalGeometries.Groups: GeometryVectorLike

"""
    sample_points_inside(geometry; n_samples=200, sampling_method=:pseudo_volume, oversampling=1.2)

Sample points from a mesh.

# Arguments
- `geometry`: A mesh geometry (`Mesh`, `FixedGeometry`, or any geometry containing mesh obstructions)
- `n_samples`: Number of points to sample (default: 200)
- `sampling_method`: `:pseudo_volume` or `:surface` (default: `:pseudo_volume`)
- `oversampling`: Oversampling factor for the pseudo-volume sampling method (default: 1.2)

# Returns
A `Vector{SVector{3, Float64}}` of sampled 3D positions.

# Sampling methods
- **`:surface`**: Blue-noise (Poisson disk) sampling on the mesh surface using sample elimination (Yuksel 2015).
- **`:pseudo_volume`**: Homogeneous surface sampling followed by a randomized inward shift along the surface normal. Points that end up outside the mesh after shifting are discarded.
"""
function sample_points_inside end

function sample_points_inside(geometry; n_samples::Int=200, sampling_method::Symbol=:pseudo_volume, oversampling::Float64=1.2)
    if n_samples <= 0
        return SVector{3, Float64}[]
    end
    fg = geometry isa FixedGeometry ? geometry : Geometries.fix(geometry)
    sample_points_inside(fg; n_samples, sampling_method, oversampling)
end

function sample_points_inside(fg::FixedGeometry; n_samples::Int=200, sampling_method::Symbol=:pseudo_volume, oversampling::Float64=1.2)
    mesh_data_list = geometry_mesh(fg)
    
    if isnothing(mesh_data_list) || isempty(mesh_data_list)
        error("No mesh found in geometry. sample_points_inside only works with mesh geometries.")
    end

    # Aggregate all triangles from all meshes in the list
    n_faces = sum(length(m.triangles) for m in mesh_data_list)
    if n_faces == 0
        error("No triangles found in mesh.")
    end
    
    face_normals = Vector{SVector{3, Float64}}(undef, n_faces)
    face_areas = Vector{Float64}(undef, n_faces)
    face_vertices = Vector{Tuple{SVector{3, Float64}, SVector{3, Float64}, SVector{3, Float64}}}(undef, n_faces)
    
    total_area = 0.0
    i = 1
    for m in mesh_data_list
        vertices = m.vertices
        for tri_idx in m.triangles
            v0, v1, v2 = vertices[tri_idx[1]], vertices[tri_idx[2]], vertices[tri_idx[3]]
            face_vertices[i] = (v0, v1, v2)
            
            edge1 = v1 - v0
            edge2 = v2 - v0
            n_unnorm = edge1 × edge2
            area = norm(n_unnorm) / 2
            if area > 0
                face_normals[i] = n_unnorm / norm(n_unnorm)
            else
                face_normals[i] = SVector(0.0, 0.0, 0.0)
            end
            face_areas[i] = area
            total_area += area
            i += 1
        end
    end

    if total_area <= 0
        error("Mesh has zero total surface area. Cannot sample points.")
    end

    max_shift = size_scale(fg)
    if !isfinite(max_shift) || max_shift <= 0
        bb = Geometries.BoundingBox(fg)
        bb_size = upper(bb) - lower(bb)
        max_shift = minimum(bb_size) / 2
    end

    n_oversample = if sampling_method == :pseudo_volume
        Int(ceil(n_samples * oversampling * 10))
    else
        n_samples * 10
    end

    pts, face_idx = _sample_surface(face_vertices, face_areas, total_area, n_oversample)

    n_target = if sampling_method == :pseudo_volume
        Int(ceil(oversampling * n_samples))
    else
        n_samples
    end
    pts, face_idx = _sample_elimination(pts, face_idx, n_target)

    if sampling_method == :pseudo_volume
        pts = _shift_inward!(pts, face_idx, face_normals, max_shift, fg)

        inside_mask = [length(isinside(fg, p)) > 0 for p in pts]
        pts = pts[inside_mask]

        if length(pts) < n_samples
            @warn "Could only sample $(length(pts)) points inside the mesh (requested $n_samples). Try rerunning with higher oversampling factor."
        end
    end

    Random.shuffle!(pts)
    pts = pts[1:min(n_samples, length(pts))]

    return pts
end

function _sample_surface(face_vertices, face_areas, total_area, n_samples)
    pts = SVector{3, Float64}[]
    face_idx = Int[]

    for (i, vertices) in enumerate(face_vertices)
        lambda = face_areas[i] / total_area * n_samples
        n_on_face = rand(Poisson(lambda))
        n_on_face == 0 && continue

        v0, v1, v2 = vertices
        e1 = v1 - v0
        e2 = v2 - v0

        for _ in 1:n_on_face
            u1, u2 = rand(2)
            if u1 + u2 > 1
                u1 = 1 - u1
                u2 = 1 - u2
            end
            push!(pts, v0 + u1 * e1 + u2 * e2)
            push!(face_idx, i)
        end
    end

    return pts, face_idx
end

function _sample_elimination(pts::Vector{SVector{3, Float64}}, face_idx::Vector{Int}, n_target::Int)
    n = length(pts)
    n <= n_target && return pts, face_idx

    mask = ones(Bool, n)
    n_to_eliminate = n - n_target
    batch_size = max(1, n_to_eliminate ÷ 50)

    while n_to_eliminate > 0
        active = findall(mask)
        active_pts = pts[active]

        tree = KDTree(active_pts)

        min_dists = Vector{Float64}(undef, length(active))
        for (j, p) in enumerate(active_pts)
            idxs, dists = knn(tree, p, 2)
            min_dists[j] = dists[2]
        end

        num_this_batch = min(batch_size, n_to_eliminate)
        remove_order = partialsortperm(min_dists, 1:num_this_batch)
        mask[active[remove_order]] .= false
        n_to_eliminate -= num_this_batch
    end

    return pts[mask], face_idx[mask]
end

function _shift_inward!(pts, face_idx, face_normals, max_shift, fg)
    n = length(pts)
    shifted = copy(pts)
    shift_amounts = rand(n) .* max_shift

    for i in 1:n
        shifted[i] = pts[i] - face_normals[face_idx[i]] * shift_amounts[i]

        if length(isinside(fg, shifted[i])) == 0
            for scale in [0.5, 0.25, 0.1, 0.05]
                shifted[i] = pts[i] - face_normals[face_idx[i]] * (shift_amounts[i] * scale)
                if length(isinside(fg, shifted[i])) > 0
                    break
                end
            end
            if length(isinside(fg, shifted[i])) == 0
                shifted[i] = pts[i]
            end
        end
    end

    return shifted
end

end
