struct Round{N} <: BaseObstruction{N}
    radius::Float64
end

has_inside(::Type{<:Round}) = true
has_single_inside(::Type{<:Round}) = true

function inside_indices(
    round::Round{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    !Base.isempty(intersection) && return intersection.inside ? [ObstructionIndex()] : ObstructionIndex[]
    isinside(round, position) ? [ObstructionIndex()] : ObstructionIndex[]
end

const InfiniteCylinder = Round{2}
const Sphere = Round{3}

InternalBoundingBox(round::Round{N}) where {N} = InternalBoundingBox{N}(round.radius)

isinside(round::Round{N}, position::SVector{N, Float64}) where {N} =
    sum(position .* position) < round.radius^2

size_scale(round::Round) = round.radius

function detect_intersection(
    round::Round{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    previous = !Base.isempty(previous_hit)
    inside = previous ? previous_hit.inside : isinside(round, start)
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

function _geometry_mesh(sphere::Sphere; nsamples=100, kwargs...)
    p = 1 / ((1 + sqrt(5)) / 2)
    vertices = [SVector{3, Float64}(v / norm(v)) for v in [[0,p,-1],[p,1,0],[-p,1,0],[0,p,1],
        [0,-p,1],[-1,0,p],[0,-p,-1],[1,0,-p],[1,0,p],[-1,0,-p],[p,-1,0],[-p,-1,0]]]
    triangles = [SVector{3, Int}(t) for t in [[3,2,1],[2,3,4],[6,5,4],[5,9,4],[8,7,1],
        [7,10,1],[12,11,5],[11,12,7],[10,6,3],[6,10,12],[9,8,2],[8,9,11],[3,6,4],[9,2,4],
        [10,3,1],[2,8,1],[12,10,7],[8,11,7],[6,12,5],[11,9,5]]]
    [_mesh_result(sphere.radius .* vertices, triangles)]
end
