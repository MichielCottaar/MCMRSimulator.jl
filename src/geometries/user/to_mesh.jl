"""
Methods that extend the `Mesh` constructor to convert different types into meshes.
"""
module ToMesh
import StaticArrays: SVector, MVector
import BSplineKit: BSplineOrder, interpolate, Derivative
import LinearAlgebra: cross, norm, ⋅
import Statistics: mean
import ..Obstructions: Mesh, value_as_vector, BendyCylinder

"""
    Mesh(bendy_cylinder)

Converts a [`BendyCylinder`](@ref) object into 
"""
function Mesh(bendy_cylinder::BendyCylinder)
    control_points = value_as_vector(bendy_cylinder.control_point)
    radii = value_as_vector(bendy_cylinder.radius)
    triangle_size = 2π * mean(radii) / bendy_cylinder.nsamples.value

    nadd = isnothing(bendy_cylinder.closed.value) ? 0 : div(bendy_cylinder.spline_order.value, length(bendy_cylinder.control_point), RoundUp)
    if iszero(nadd)
        full_control_points = SVector{3}.(control_points)
        full_radii = Float64.(radii)
    else
        @assert !isnothing(bendy_cylinder.repeats.value)
        full_control_points = [SVector{3}(cp .+ bendy_cylinder.repeats.value .* bendy_cylinder.closed.value .* index_add) for index_add in -nadd:nadd for cp in control_points]
        full_radii = [r for _ in -nadd:nadd for r in radii]
    end

    dist = [0, cumsum(norm.(full_control_points[2:end] .- full_control_points[1:end-1]))...]
    dist .-= dist[nadd * length(control_points) + 1]  # set first of original control points to reference point with distance of 0

    fpos = interpolate(dist, full_control_points, BSplineOrder(bendy_cylinder.spline_order.value))
    frad = interpolate(dist, full_radii, BSplineOrder(bendy_cylinder.spline_order.value))

    nsamples = bendy_cylinder.nsamples.value
    theta = range(start=0, stop=2π, length=nsamples + 1)[1:end-1]
    if iszero(nadd)
        dist_control = dist
    else
        dist_control = dist[nadd * length(control_points) + 1: (nadd + 1) * length(control_points) + 1]
    end
    dist_eval = [0., [
        d
        for (d1, d2) in zip(dist_control[1:end-1], dist_control[2:end])
        for d in range(start=d1, stop=d2, length=max(Int(div(d2 - d1, triangle_size, RoundUp)), 2))[2:end]
    ]...]

    vertices = SVector{3, Float64}[]
    triangles = SVector{3, Int}[]

    fdirection = Derivative(1) * fpos

    # start off vec_theta0 in arbitrary direction
    vec_theta0_try = cross(fdirection(0), SVector{3}((1., 0., 0.)))    
    if all(iszero.(vec_theta0_try))
        vec_theta0_try = cross(fdirection(0), SVector{3}((0., 1., 0.)))    
    end
    vec_theta0 = MVector{3}(vec_theta0_try ./ norm(vec_theta0_try))

    function get_vertices(dist)
        # adjust vec_theta0 to be orthogonal to current direction
        direction_unnorm = fdirection(dist)
        direction = direction_unnorm ./ norm(direction_unnorm)
        vec_theta0 .-= direction .* (direction ⋅ vec_theta0)
        vec_theta0 ./= norm(vec_theta0)
        centroid = fpos(dist)
        radius_uncorrected = frad(dist)
        # correct radius to preserve final volume
        half_theta_step = π / nsamples
        radius = radius_uncorrected * (half_theta_step / (sin(half_theta_step) * cos(half_theta_step)))
        last_vec = cross(direction, vec_theta0)
        return [SVector{3}((radius * sin(t)) .* last_vec .+ (radius * cos(t)) .* vec_theta0 .+ centroid) for t in theta]
    end

    vertices = vcat(get_vertices.(dist_eval)...)

    single_triangle_set = [
        SVector{3}(nsamples, 1, nsamples * 2),
        [SVector{3}((add, 1 + add, nsamples + add)) for add in 1:nsamples-1]...,
        SVector{3}(1, nsamples * 2, nsamples + 1),
        [SVector{3}((add + 1, nsamples + add, nsamples + add + 1)) for add in 1:nsamples-1]...,
    ]

    triangles = [t .+ index * nsamples for index in 0:length(dist_eval)-2 for t in single_triangle_set]

    mesh_kwargs = Dict{Symbol, Any}(
        :vertices => vertices,
        :triangles => triangles,
        :number => length(triangles),
    )
    for symbol in Mesh(number=0).unique_keys
        if symbol in keys(mesh_kwargs)
            continue
        end
        mesh_kwargs[symbol] = getproperty(bendy_cylinder, symbol).value
    end

    return Mesh(; mesh_kwargs...)
end

end