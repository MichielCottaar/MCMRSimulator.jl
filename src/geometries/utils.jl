module Utils

import LinearAlgebra: cross, ⋅
import Quickhull: facets, quickhull
import StaticArrays: SVector

function volume_conserving_cylinder_radius(radius, nsamples)
    half_theta_step = π / nsamples
    radius * sqrt(half_theta_step / (sin(half_theta_step) * cos(half_theta_step)))
end

function volume_conserving_sphere_radius(radius, vertices, triangles)
    polyhedron_volume = abs(sum(
        vertices[triangle[1]] ⋅ cross(vertices[triangle[2]], vertices[triangle[3]]) / 6
        for triangle in triangles
    ))
    radius * (4π / (3 * polyhedron_volume))^(1 / 3)
end

"""
    sphere_mesh(nsamples)

Generate an approximately uniform triangular mesh on the unit sphere using a
Fibonacci lattice and a three-dimensional convex hull.
"""
function sphere_mesh(nsamples::Int)
    nsamples >= 4 || throw(ArgumentError("nsamples must be at least 4"))

    golden_angle = π * (3 - sqrt(5))
    vertices = [
        let
            z = 1 - 2 * (index - 0.5) / nsamples
            radial_distance = sqrt(1 - z^2)
            theta = golden_angle * (index - 1)
            SVector(radial_distance * cos(theta), radial_distance * sin(theta), z)
        end
        for index in 1:nsamples
    ]

    hull = quickhull(vertices)
    triangles = [SVector{3, Int}(Tuple(face)) for face in facets(hull)]
    triangles = [
        let
            triangle_normal = cross(vertices[triangle[2]] - vertices[triangle[1]],
                vertices[triangle[3]] - vertices[triangle[1]])
            dot_product = (sum(vertices[index] for index in triangle) / 3) ⋅ triangle_normal
            dot_product < 0 ? SVector{3, Int}(triangle[1], triangle[3], triangle[2]) : triangle
        end
        for triangle in triangles
    ]

    return vertices, triangles
end

end
