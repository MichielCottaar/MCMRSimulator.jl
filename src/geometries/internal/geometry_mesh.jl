"""Render-mesh extraction for fixed physical geometries."""
module GeometryMeshes

import StaticArrays: SVector
import LinearAlgebra: cross, norm
import ..Internal: FixedGeometry, geometry_mesh
import ...BoundingBoxes: BoundingBox
import ..InternalBoundingBoxes
import ..PhysicalGeometries: Transparent, transparent_geometry, InternalBoundingBox
import ..PhysicalGeometries.Groups: GeometryVectorLike, GeometryTuple
import ..PhysicalGeometries.Transformations: Shift, Scale, Rotate, backward, forward
import ..PhysicalGeometries.Repeats: Repeat
import ..PhysicalGeometries.Meshes: Mesh
import ..PhysicalGeometries.BaseObstructions: InfiniteWall, InfiniteCylinder, Sphere

function _mesh_result(vertices, triangles)
    (vertices=vertices, triangles=triangles)
end

_internal_box(::Nothing) = nothing
_internal_box(box::BoundingBox) = InternalBoundingBox(box)
_internal_box(box::InternalBoundingBox) = box

function geometry_mesh(
    geometry;
    height=nothing,
    nsamples=100,
    bounding_box=nothing,
)
    _geometry_mesh(
        geometry;
        height,
        nsamples,
        bounding_box=_internal_box(bounding_box),
    )
end

geometry_mesh(geometry, bounding_box; kwargs...) =
    geometry_mesh(geometry; bounding_box, kwargs...)

function _geometry_mesh(mesh::Mesh; kwargs...)
    number_of_triangles = mesh.first_index_of_gap - 1
    [_mesh_result(mesh.vertices, mesh.indices[1:number_of_triangles])]
end

_geometry_mesh(::InfiniteWall; kwargs...) = [0.]

function _geometry_mesh(cylinder::InfiniteCylinder; nsamples=100, kwargs...)
    nsamples >= 3 || throw(ArgumentError("nsamples must be at least 3"))
    [[
        SVector{2, Float64}(
            cylinder.radius * cos(2π * index / nsamples),
            cylinder.radius * sin(2π * index / nsamples),
        )
        for index in 0:(nsamples - 1)
    ]]
end

function _icosahedron()
    golden_ratio = (1 + sqrt(5)) / 2
    p = 1 / golden_ratio
    vertices = [v / norm(v) for v in [
        [0, p, -1], [p, 1, 0], [-p, 1, 0], [0, p, 1],
        [0, -p, 1], [-1, 0, p], [0, -p, -1], [1, 0, -p],
        [1, 0, p], [-1, 0, -p], [p, -1, 0], [-p, -1, 0],
    ]]
    triangles = [
        [3, 2, 1], [2, 3, 4], [6, 5, 4], [5, 9, 4], [8, 7, 1],
        [7, 10, 1], [12, 11, 5], [11, 12, 7], [10, 6, 3], [6, 10, 12],
        [9, 8, 2], [8, 9, 11], [3, 6, 4], [9, 2, 4], [10, 3, 1],
        [2, 8, 1], [12, 10, 7], [8, 11, 7], [6, 12, 5], [11, 9, 5],
    ]
    SVector{3, Float64}.(vertices), SVector{3, Int}.(triangles)
end

function _icosahedron(subdivisions::Int)
    vertices, base_triangles = _icosahedron()
    isone(subdivisions) && return vertices, base_triangles

    edge_map = Dict{Tuple{Int, Int}, Vector{Int}}()
    for triangle in base_triangles
        for (i1, i2) in ((triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1]))
            haskey(edge_map, (i1, i2)) && continue
            edge_map[(i1, i2)] = length(vertices) .+ (1:(subdivisions - 1))
            edge_map[(i2, i1)] = reverse(edge_map[(i1, i2)])
            append!(vertices, [
                vertices[i1] .* (1 - distance) .+ vertices[i2] .* distance
                for distance in (0:1 / (subdivisions + 1):1)[2:end-1]
            ])
        end
    end

    triangles = SVector{3, Int}[]
    for triangle in base_triangles
        vertex_indices = Dict{Tuple{Int, Int}, Int}()
        for layer in 1:(subdivisions + 1)
            for index in 1:(subdivisions + 2 - layer)
                if layer == 1
                    vertex_indices[(layer, index)] = if index == 1
                        triangle[1]
                    elseif index == subdivisions + 1
                        triangle[2]
                    else
                        edge_map[(triangle[1], triangle[2])][index - 1]
                    end
                elseif layer == subdivisions + 1
                    vertex_indices[(layer, index)] = triangle[3]
                elseif index == 1
                    vertex_indices[(layer, index)] = edge_map[(triangle[1], triangle[3])][layer - 1]
                elseif index == subdivisions + 2 - layer
                    vertex_indices[(layer, index)] = edge_map[(triangle[2], triangle[3])][layer - 1]
                else
                    v1 = vertices[edge_map[(triangle[1], triangle[3])][layer - 1]]
                    v2 = vertices[edge_map[(triangle[2], triangle[3])][layer - 1]]
                    distance = (index - 1) / (subdivisions + 1 - layer)
                    push!(vertices, (1 - distance) .* v1 .+ distance .* v2)
                    vertex_indices[(layer, index)] = length(vertices)
                end
            end
        end
        for layer in 1:subdivisions
            for index in 1:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}(
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer, index + 1)],
                    vertex_indices[(layer + 1, index)],
                ))
            end
        end
        for layer in 1:(subdivisions - 1)
            for index in 2:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}(
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer + 1, index)],
                    vertex_indices[(layer + 1, index - 1)],
                ))
            end
        end
    end
    [v / norm(v) for v in vertices], triangles
end

function _geometry_mesh(sphere::Sphere; nsamples=100, kwargs...)
    subdivisions = Int(ceil(sqrt(nsamples / 20)))
    vertices, triangles = _icosahedron(subdivisions)
    vertices = [sphere.radius .* vertex for vertex in vertices]
    [_mesh_result(vertices, triangles)]
end

function _transform_mesh(transformation, mesh_data)
    [
        _mesh_result(
            [backward(transformation, vertex) for vertex in mesh.vertices],
            mesh.triangles,
        )
        for mesh in mesh_data
    ]
end

function _translate_native(value, displacement)
    value isa Real && return value + displacement[1]
    value isa SVector && return value + displacement
    value isa NamedTuple && return _mesh_result(value.vertices .+ Ref(displacement), value.triangles)
    [_translate_native(item, displacement) for item in value]
end

function _transform_native(transformation, value)
    value isa Real && return backward(transformation, SVector{1, Float64}(value))[1]
    value isa SVector && return backward(transformation, value)
    value isa NamedTuple && return _transform_mesh(transformation, [value])[1]
    [_transform_native(transformation, item) for item in value]
end

function _height(height, bounding_box, axis)
    !isnothing(height) && return Float64(height)
    isnothing(bounding_box) && return 1.
    2 * sum(abs.(axis) .* InternalBoundingBoxes.half_size(bounding_box))
end

function _wall_basis(normal)
    reference = abs(normal[1]) < 0.9 ? SVector(1., 0., 0.) : SVector(0., 1., 0.)
    first = cross(normal, reference)
    first ./ norm(first), cross(normal, first ./ norm(first))
end

function _deproject_walls(positions, transformation::Rotate{3, 1}; height, bounding_box)
    normal = SVector{3, Float64}(transformation.matrix[1, :])
    first, second = _wall_basis(normal)
    first_height = _height(height, bounding_box, first)
    second_height = _height(height, bounding_box, second)
    result = NamedTuple[]
    for position in positions
        center = transformation.matrix' * SVector{1, Float64}(position)
        vertices = [
            center + first_height / 2 * first + second_height / 2 * second,
            center - first_height / 2 * first + second_height / 2 * second,
            center + first_height / 2 * first - second_height / 2 * second,
            center - first_height / 2 * first - second_height / 2 * second,
        ]
        push!(result, _mesh_result(vertices, [SVector(1, 2, 3), SVector(4, 3, 2)]))
    end
    result
end

function _deproject_circles(circles, transformation::Rotate{3, 2}; height, bounding_box)
    first = SVector{3, Float64}(transformation.matrix[1, :])
    second = SVector{3, Float64}(transformation.matrix[2, :])
    axis = cross(first, second)
    axis = axis ./ norm(axis)
    extrusion = _height(height, bounding_box, axis)
    result = NamedTuple[]
    for circle in circles
        points = [transformation.matrix' * point for point in circle]
        vertices = vcat(points .+ Ref(axis * extrusion / 2), points .- Ref(axis * extrusion / 2))
        count = length(points)
        triangles = SVector{3, Int}[]
        for index in 1:count
            next = mod1(index + 1, count)
            push!(triangles, SVector(index, next, count + index))
            push!(triangles, SVector(next, count + next, count + index))
        end
        push!(result, _mesh_result(vertices, triangles))
    end
    result
end

function _geometry_mesh_preserving(transformation, geometry; bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    _transform_native(
        transformation,
        _geometry_mesh(geometry; kwargs..., bounding_box=local_box),
    )
end

_geometry_mesh(transformation::Shift, geometry; kwargs...) =
    _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Scale, geometry; kwargs...) =
    _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Rotate{N, N}, geometry; kwargs...) where N =
    _geometry_mesh_preserving(transformation, geometry; kwargs...)

_geometry_mesh(transformation::Shift; kwargs...) =
    _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Scale; kwargs...) =
    _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Rotate; kwargs...) =
    _geometry_mesh(transformation, transformation.geometry; kwargs...)

function _geometry_mesh(transformation::Rotate{3, 1}, geometry; height=nothing, nsamples=100, bounding_box=nothing)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    positions = _geometry_mesh(
        geometry;
        height,
        nsamples,
        bounding_box=local_box,
    )
    _deproject_walls(positions, transformation; height, bounding_box)
end

function _geometry_mesh(transformation::Rotate{3, 2}, geometry; height=nothing, nsamples=100, bounding_box=nothing)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    circles = _geometry_mesh(
        geometry;
        height,
        nsamples,
        bounding_box=local_box,
    )
    _deproject_circles(circles, transformation; height, bounding_box)
end

function _geometry_mesh(geometry::GeometryVectorLike; kwargs...)
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])
end

function _geometry_mesh(geometry::GeometryTuple; kwargs...)
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])
end

function _geometry_mesh(repeat::Repeat{N}; bounding_box=nothing, kwargs...) where N
    child_box = InternalBoundingBox(repeat.geometry)
    lower_index, upper_index = if isnothing(bounding_box)
        zeros(Int, N), zeros(Int, N)
    else
        (
            floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats),
            ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats),
        )
    end
    child = _geometry_mesh(repeat.geometry; kwargs...)
    result = Any[]
    for index in Iterators.product((lower_index[i]:upper_index[i] for i in 1:N)...)
        displacement = SVector{N, Float64}(index) .* repeat.repeats
        if !isnothing(bounding_box)
            shifted_lower = InternalBoundingBoxes.lower(child_box) + displacement
            shifted_upper = InternalBoundingBoxes.upper(child_box) + displacement
            all(shifted_lower .<= InternalBoundingBoxes.upper(bounding_box)) || continue
            all(shifted_upper .>= InternalBoundingBoxes.lower(bounding_box)) || continue
        end
        append!(result, [_translate_native(value, displacement) for value in child])
    end
    result
end

_geometry_mesh(wrapper::Transparent; kwargs...) =
    _geometry_mesh(transparent_geometry(wrapper); kwargs...)
_geometry_mesh(geometry::FixedGeometry; kwargs...) =
    _geometry_mesh(geometry.geometry; kwargs...)

end
