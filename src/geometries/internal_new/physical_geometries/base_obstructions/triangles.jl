import LinearAlgebra: cross, norm, ⋅

struct FullTriangle <: BaseObstruction{3}
    a::SVector{3, Float64}
    b::SVector{3, Float64}
    c::SVector{3, Float64}
end

has_inside(::Type{FullTriangle}) = false
inside_indices(::FullTriangle, ::SVector{3, Float64}) = ObstructionIndex[]

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

normal(triangle::FullTriangle) = normal(triangle.a, triangle.b, triangle.c)
normal(a::AbstractVector, b::AbstractVector, c::AbstractVector) = begin
    unnormalized = cross(b - a, c - a)
    unnormalized ./ norm(unnormalized)
end

function detect_intersection(
    triangle::FullTriangle,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
)
    !Base.isempty(previous_hit) && return Intersection{3}()
    triangle_normal = normal(triangle)
    plane_distance = triangle_normal ⋅ triangle.a
    start_distance = triangle_normal ⋅ start
    destination_distance = triangle_normal ⋅ destination
    abs(start_distance - destination_distance) < 1e-8 && return Intersection{3}()

    distance = (plane_distance - start_distance) / (destination_distance - start_distance)
    (distance < 0 || distance > 1) && return Intersection{3}()
    intersection_point = distance .* destination .+ (1 - distance) .* start

    for (dimension, next_dimension) in ((1, 2), (2, 3), (3, 1))
        edge = triangle[next_dimension] - triangle[dimension]
        to_point = intersection_point - triangle[dimension]
        cross(edge, to_point) ⋅ triangle_normal < 0 && return Intersection{3}()
    end

    inside = destination_distance > start_distance
    return Intersection(distance, inside ? -triangle_normal : triangle_normal, inside, ObstructionIndex(), false)
end
