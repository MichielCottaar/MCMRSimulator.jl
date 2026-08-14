"""Surface sampling for the physical-geometry hierarchy."""
module SurfaceSamplings

import StaticArrays: SVector
import Distributions: Poisson
import LinearAlgebra: cross, norm, nullspace
import Random: rand

import ...InternalBoundingBoxes
import ...InternalBoundingBoxes: InternalBoundingBox
import ...Properties: GeometryProperties, GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties
import ..PhysicalGeometries: PhysicalGeometry, surface_sampling
import ..Groups: GeometryVectorLike, GeometryTuple
import ..Transformations: Transformation, Shift, Scale, Rotate, forward, backward, backward_normal
import ..Repeats: Repeat
import ..Transparents: Transparent, transparent_geometry
import ..BaseObstructions: InfiniteWall, Round, OverlappingRound, FullTriangle, normal
import ..Meshes
import ..Meshes: Mesh

function _combine(draws, ::Val{N}) where {N}
    draws = collect(draws)
    positions = reduce(vcat, (draw[1] for draw in draws); init=SVector{N, Float64}[])
    normals = reduce(vcat, (draw[2] for draw in draws); init=SVector{N, Float64}[])
    positions, normals
end

_density_child(density::GeometryLeafProperties, ::Int) = density
_density_child(density::GeometryVectorProperties, index::Int) = density.properties[index]
_density_child(density::GeometryTupleProperties, index::Int) = density.properties[index]

function surface_sampling(
    geometry::GeometryVectorLike{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    _combine(
        (
            surface_sampling(
                child,
                _density_child(density, index),
                bounding_box,
                scale_density,
            )
            for (index, child) in enumerate(geometry)
        ),
        Val(N),
    )
end

function surface_sampling(
    geometry::GeometryTuple{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    _combine(
        (
            surface_sampling(
                child,
                _density_child(density, index),
                bounding_box,
                scale_density,
            )
            for (index, child) in enumerate(geometry)
        ),
        Val(N),
    )
end

function surface_sampling(
    transformation::Shift{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    positions, normals = surface_sampling(
        transformation.geometry,
        density,
        backward(transformation, bounding_box),
        scale_density,
    )
    positions = [backward(transformation, position) for position in positions]
    positions, normals
end

function surface_sampling(
    transformation::Scale{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    positions, normals = surface_sampling(
        transformation.geometry,
        density,
        backward(transformation, bounding_box),
        scale_density * transformation.scale^(N - 1),
    )
    positions = [backward(transformation, position) for position in positions]
    positions, normals
end

function _deproject_positions(
    transformation::Rotate{N, M},
    positions,
    bounding_box::InternalBoundingBox{N},
) where {N, M}
    basis = nullspace(Matrix(transformation.matrix))
    center = InternalBoundingBoxes.center(bounding_box)
    half_size = InternalBoundingBoxes.half_size(bounding_box)
    null_center = basis' * center
    null_half_size = abs.(basis)' * half_size
    [
        transformation.matrix' * position + basis * (
            null_center + (2 .* rand(length(null_center)) .- 1) .* null_half_size
        )
        for position in positions
    ]
end

function _projected_scale(
    transformation::Rotate{N, M},
    bounding_box::InternalBoundingBox{N},
) where {N, M}
    basis = nullspace(Matrix(transformation.matrix))
    half_size = InternalBoundingBoxes.half_size(bounding_box)
    prod(2 .* (abs.(basis)' * half_size))
end

function surface_sampling(
    transformation::Rotate{N, M},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N, M}
    if N == M
        positions, normals = surface_sampling(
            transformation.geometry,
            density,
            backward(transformation, bounding_box),
            scale_density,
        )
        return [backward(transformation, position) for position in positions],
            [backward_normal(transformation, normal) for normal in normals]
    end

    positions, normals = surface_sampling(
        transformation.geometry,
        density,
        forward(transformation, bounding_box),
        scale_density * _projected_scale(transformation, bounding_box),
    )
    positions = _deproject_positions(transformation, positions, bounding_box)
    normals = [backward_normal(transformation, normal) for normal in normals]
    keep = [
        all(position .>= InternalBoundingBoxes.lower(bounding_box)) &&
        all(position .<= InternalBoundingBoxes.upper(bounding_box))
        for position in positions
    ]
    positions[keep], normals[keep]
end

function surface_sampling(
    wrapper::Transparent{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    surface_sampling(transparent_geometry(wrapper), density, bounding_box, scale_density)
end

function surface_sampling(
    repeat::Repeat{N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    child_box = InternalBoundingBox(repeat.geometry)
    lower = floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats)
    upper = ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats)
    draws = (
        let displacement = SVector{N, Float64}(indices) .* repeat.repeats
            positions, normals = surface_sampling(
                repeat.geometry,
                density,
                InternalBoundingBoxes.shift(bounding_box, -displacement),
                scale_density,
            )
            ([position + displacement for position in positions], normals)
        end
        for indices in Iterators.product((lower[index]:upper[index] for index in 1:N)...)
    )
    positions, normals = _combine(draws, Val(N))
    keep = [
        all(position .>= InternalBoundingBoxes.lower(bounding_box)) &&
        all(position .<= InternalBoundingBoxes.upper(bounding_box))
        for position in positions
    ]
    positions[keep], normals[keep]
end

function _filter_to_box(positions, normals, bounding_box)
    keep = [
        all(position .>= InternalBoundingBoxes.lower(bounding_box)) &&
        all(position .<= InternalBoundingBoxes.upper(bounding_box))
        for position in positions
    ]
    positions[keep], normals[keep]
end

function _round_sampling(round::Round{N}, density, scale_density) where {N}
    effective_density = density * scale_density
    if N == 2
        surface = 2π * round.radius
        nspins = rand(Poisson(surface * effective_density))
        theta = rand(nspins) .* 2π
        normals = [SVector{2, Float64}(cos(t), sin(t)) for t in theta]
        return normals .* (-round.radius), normals
    end

    surface = 4π * round.radius^2
    nspins = rand(Poisson(surface * effective_density))
    normals = [
        let z = rand(Float64) * 2 - 1, r = sqrt(1 - z^2), theta = rand(Float64) * 2π
            s, c = sincos(theta)
            SVector{3, Float64}(r * s, r * c, z)
        end
        for _ in 1:nspins
    ]
    normals .* (-round.radius), normals
end

function surface_sampling(
    round::Round{N},
    density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    _filter_to_box(_round_sampling(round, density.value, scale_density)..., bounding_box)
end

function surface_sampling(
    round::OverlappingRound{N},
    density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density,
) where {N}
    surface_sampling(Round{N}(round.radius), density, bounding_box, scale_density)
end

function surface_sampling(
    ::InfiniteWall,
    density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{1},
    scale_density,
)
    nspins = rand(Poisson(density.value * scale_density))
    _filter_to_box(
        fill(zero(SVector{1, Float64}), nspins),
        fill(SVector{1, Float64}(1.0), nspins),
        bounding_box,
    )
end

function _triangle_sampling(triangle::FullTriangle, density, scale_density)
    edge_1 = triangle.b - triangle.a
    edge_2 = triangle.c - triangle.a
    surface = norm(cross(edge_1, edge_2)) / 2
    nspins = rand(Poisson(surface * density * scale_density))
    draw_position() = begin
        u1, u2 = rand(2)
        if u1 + u2 > 1
            u1 = 1 - u1
            u2 = 1 - u2
        end
        SVector{3, Float64}(triangle.a + u1 * edge_1 + u2 * edge_2)
    end
    positions = [draw_position() for _ in 1:nspins]
    positions, fill(-normal(triangle), nspins)
end

function surface_sampling(
    triangle::FullTriangle,
    density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{3},
    scale_density,
)
    _filter_to_box(_triangle_sampling(triangle, density.value, scale_density)..., bounding_box)
end

function surface_sampling(
    mesh::Mesh,
    density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{3},
    scale_density,
)
    draws = (_triangle_sampling(Meshes.triangle(mesh, index), density.value, scale_density) for index in eachindex(mesh.indices))
    _filter_to_box(_combine(draws, Val(3))..., bounding_box)
end

end
