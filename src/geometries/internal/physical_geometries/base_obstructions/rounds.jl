import ....Utils.ToMesh: icosahedron

struct Round{N} <: BaseObstruction{N}
    radius::Float64
end

has_inside(::Type{<:Round}) = true
has_single_inside(::Type{<:Round}) = true

function isinside_single(
    round::Round{N},
    position::SVector{N, Float64},
    intersection::Intersection{3}=Intersection{3}(),
) where {N}
    !Base.isempty(intersection) && return intersection.inside
    return sum(position .* position) < round.radius^2
end

const InfiniteCylinder = Round{2}
const Sphere = Round{3}

InternalBoundingBox(round::Round{N}) where {N} = InternalBoundingBox{N}(round.radius)

size_scale(round::Round) = round.radius

function detect_intersection(
    round::Round{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    previous = !Base.isempty(previous_hit)
    inside = previous ? previous_hit.inside : isinside_single(round, start)
    !inside && previous && return Intersection{N}()
    difference = destination - start
    a = sum(difference .* difference)
    b = sum(2 .* start .* difference)
    c = sum(start .* start)
    determinant = b * b - 4 * a * (c - round.radius^2)
    determinant < 0 && return Intersection{N}()

    solution = (inside ? -b + sqrt(determinant) : -b - sqrt(determinant)) / (2 * a)
    (solution <= 0 || solution > 1) && return Intersection{N}()

    normal = (solution .* destination .+ (1 - solution) .* start) ./ round.radius
    return Intersection(solution, inside ? -normal : normal, inside, ObstructionIndex(), false)
end

function _round_sampling(round::Round{N}, density, scale_density) where {N}
    effective_density = density * scale_density
    if N == 2
        nspins = rand(Poisson(2π * round.radius * effective_density))
        theta = rand(nspins) .* 2π
        normals = [SVector{2, Float64}(cos(t), sin(t)) for t in theta]
        return normals .* (-round.radius), normals
    end
    nspins = rand(Poisson(4π * round.radius^2 * effective_density))
    normals = [let z = rand(Float64) * 2 - 1, r = sqrt(1 - z^2), theta = rand(Float64) * 2π
        s, c = sincos(theta)
        SVector{3, Float64}(r * s, r * c, z)
    end for _ in 1:nspins]
    normals .* (-round.radius), normals
end

function surface_sampling(
    round::Round{N}, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{N}, scale_density,
) where {N}
    _round_sampling(round, density.value, scale_density)
end

function _geometry_mesh(cylinder::InfiniteCylinder; nsamples=100, kwargs...)
    nsamples >= 3 || throw(ArgumentError("nsamples must be at least 3"))
    [[SVector{2, Float64}(cylinder.radius * cos(2π * index / nsamples),
        cylinder.radius * sin(2π * index / nsamples)) for index in 0:(nsamples - 1)]]
end

function _geometry_mesh(sphere::Sphere; nsamples, kwargs...)
    subdivisions = Int(ceil(sqrt(nsamples / 20)))
    base_vertices, triangles = icosahedron(subdivisions)
    [_mesh_result(sphere.radius .* base_vertices, triangles)]
end
