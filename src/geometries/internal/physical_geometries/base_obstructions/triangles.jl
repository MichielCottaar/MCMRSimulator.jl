import LinearAlgebra: cross, norm, ⋅

struct FullTriangle <: BaseObstruction{3}
    a::SVector{3, Float64}
    b::SVector{3, Float64}
    c::SVector{3, Float64}
end

has_inside(::Type{FullTriangle}) = false
has_single_inside(::Type{FullTriangle}) = false

Base.getindex(triangle::FullTriangle, index::Int) = (triangle.a, triangle.b, triangle.c)[index]

function FullTriangle(a::AbstractVector, b::AbstractVector, c::AbstractVector)
    FullTriangle(SVector{3, Float64}(a), SVector{3, Float64}(b), SVector{3, Float64}(c))
end

function InternalBoundingBox(triangle::FullTriangle)
    lower = SVector(
        min(triangle.a[1], triangle.b[1], triangle.c[1]),
        min(triangle.a[2], triangle.b[2], triangle.c[2]),
        min(triangle.a[3], triangle.b[3], triangle.c[3]),
    )
    upper = SVector(
        max(triangle.a[1], triangle.b[1], triangle.c[1]),
        max(triangle.a[2], triangle.b[2], triangle.c[2]),
        max(triangle.a[3], triangle.b[3], triangle.c[3]),
    )
    InternalBoundingBox((upper - lower) / 2, (upper + lower) / 2)
end

radius(triangle::FullTriangle) = maximum(
    InternalBoundingBoxes.upper(InternalBoundingBox(triangle)) -
    InternalBoundingBoxes.lower(InternalBoundingBox(triangle)),
) / 2

triangle_size(triangle::FullTriangle) = triangle_size(triangle.a, triangle.b, triangle.c)
triangle_size(a, b, c) = norm(cross(b - a, c - a)) / 2

normal(triangle::FullTriangle) = normal(triangle.a, triangle.b, triangle.c)
normal(a::AbstractVector, b::AbstractVector, c::AbstractVector) = begin
    unnormalized = cross(b - a, c - a)
    unnormalized ./ norm(unnormalized)
end

function find_intersection(
    triangle::FullTriangle,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit=nothing,
)
    !isnothing(previous_hit) && return nothing
    triangle_normal = normal(triangle)
    plane_distance = triangle_normal ⋅ triangle.a
    start_distance = triangle_normal ⋅ start
    destination_distance = triangle_normal ⋅ destination
    abs(start_distance - destination_distance) < 1e-8 && return nothing

    distance = (plane_distance - start_distance) / (destination_distance - start_distance)
    (distance < 0 || distance > 1) && return nothing
    intersection_point = distance .* destination .+ (1 - distance) .* start

    for (dimension, next_dimension) in ((1, 2), (2, 3), (3, 1))
        edge = triangle[next_dimension] - triangle[dimension]
        to_point = intersection_point - triangle[dimension]
        cross(edge, to_point) ⋅ triangle_normal < 0 && return nothing
    end

    inside = destination_distance > start_distance
    return (inside, distance)
end

function get_intersection_params(
    triangle::FullTriangle,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    intersection::Tuple,
)
    inside, _ = intersection
    triangle_normal = normal(triangle)
    (
        inside=inside,
        normal=inside ? -triangle_normal : triangle_normal,
        hit_gap=false,
    )
end

function surface_sampling(
    triangle::FullTriangle, density::GeometryLeafProperties,
    scale_density,
)
    edge_1 = triangle.b - triangle.a
    edge_2 = triangle.c - triangle.a
    surface = norm(cross(edge_1, edge_2)) / 2
    nspins = rand(Poisson(surface * density.value * scale_density))
    draw_position() = begin
        u1, u2 = rand(2)
        if u1 + u2 > 1
            u1 = 1 - u1
            u2 = 1 - u2
        end
        SVector{3, Float64}(triangle.a + u1 * edge_1 + u2 * edge_2)
    end
    positions = [draw_position() for _ in 1:nspins]
    positions
end
