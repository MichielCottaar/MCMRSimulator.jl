import ....Utils: sphere_mesh, volume_conserving_cylinder_radius, volume_conserving_sphere_radius
import LinearAlgebra: norm

struct Round{N} <: BaseObstruction{N}
    radius::Float64
end

distance_to_surface(round::Round{N}, position::SVector{N, Float64}) where {N} = abs(norm(position) - round.radius)

has_inside(::Type{<:Round}) = true
has_single_inside(::Type{<:Round}) = true

function isinside_single(
    round::Round{N},
    position::SVector{N, Float64},
    previous_intersection=nothing,
) where {N}
    !isnothing(previous_intersection) && return previous_intersection[1]
    return sum(position .* position) < round.radius^2
end

const InfiniteCylinder = Round{2}
const Sphere = Round{3}

InternalBoundingBox(round::Round{N}) where {N} = InternalBoundingBox{N}(round.radius)

size_scale(round::Round) = round.radius

function find_intersection(
    round::Round{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
) where {N}
    previous = !isnothing(previous_hit)
    inside = previous ? previous_hit[1] : isinside_single(round, start)
    !inside && previous && return nothing
    difference = destination - start
    a = sum(difference .* difference)
    b = sum(2 .* start .* difference)
    c = sum(start .* start)
    determinant = b * b - 4 * a * (c - round.radius^2)
    determinant < 0 && return nothing

    solution = (inside ? -b + sqrt(determinant) : -b - sqrt(determinant)) / (2 * a)
    (solution <= 0 || solution > 1) && return nothing

    return (inside, solution)
end

function get_intersection_params(
    round::Round{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    intersection::Tuple,
) where {N}
    inside, solution = intersection
    normal = (solution .* destination .+ (1 - solution) .* start) ./ round.radius
    (
        inside=inside,
        normal=inside ? -normal : normal,
        hit_gap=false,
    )
end

function _round_sampling(round::Round{N}, density, scale_density) where {N}
    effective_density = density * scale_density
    if N == 2
        nspins = rand(Poisson(2π * round.radius * effective_density))
        theta = rand(nspins) .* 2π
        normals = [SVector{2, Float64}(cos(t), sin(t)) for t in theta]
        return normals .* (-round.radius)
    end
    nspins = rand(Poisson(4π * round.radius^2 * effective_density))
    normals = [let z = rand(Float64) * 2 - 1, r = sqrt(1 - z^2), theta = rand(Float64) * 2π
        s, c = sincos(theta)
        SVector{3, Float64}(r * s, r * c, z)
    end for _ in 1:nspins]
    normals .* (-round.radius)
end

function surface_sampling(
    round::Round{N}, density::GeometryLeafProperties,
    scale_density,
) where {N}
    _round_sampling(round, density.value, scale_density)
end

function _geometry_mesh(cylinder::InfiniteCylinder; nsamples=100, kwargs...)
    nsamples >= 3 || throw(ArgumentError("nsamples must be at least 3"))
    radius = volume_conserving_cylinder_radius(cylinder.radius, nsamples)
    [[SVector{2, Float64}(radius * cos(2π * index / nsamples),
        radius * sin(2π * index / nsamples)) for index in 0:(nsamples - 1)]]
end

function _geometry_mesh(sphere::Sphere; nsamples, kwargs...)
    base_vertices, triangles = sphere_mesh(nsamples)
    radius = volume_conserving_sphere_radius(sphere.radius, base_vertices, triangles)
    [_mesh_result(radius .* base_vertices, triangles)]
end
