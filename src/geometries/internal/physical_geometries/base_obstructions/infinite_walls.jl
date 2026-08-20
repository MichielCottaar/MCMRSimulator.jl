struct InfiniteWall <: BaseObstruction{1} end

distance_to_surface(::InfiniteWall, position::SVector{1, Float64}) = abs(position[1])

has_inside(::Type{InfiniteWall}) = false
has_single_inside(::Type{InfiniteWall}) = false

const _negative_wall_normal = SVector(-1.0)
const _positive_wall_normal = SVector(1.0)

InternalBoundingBox(::InfiniteWall) = InternalBoundingBox{1}(SVector(0.0), SVector(0.0))

function find_intersection(
    ::InfiniteWall,
    start::SVector{1, Float64},
    destination::SVector{1, Float64},
    previous_hit=nothing,
)
    !isnothing(previous_hit) && return nothing
    origin = start[1]
    destination_position = destination[1]
    origin * destination_position > 0 && return nothing
    inside = origin > 0
    return (inside, abs(origin) / abs(origin - destination_position))
end

function get_intersection_params(
    ::InfiniteWall,
    start::SVector{1, Float64},
    destination::SVector{1, Float64},
    intersection::Tuple,
)
    inside, _ = intersection
    (
        inside=inside,
        normal=inside ? _positive_wall_normal : _negative_wall_normal,
        hit_gap=false,
    )
end

function surface_sampling(
    ::InfiniteWall, density::GeometryLeafProperties,
    scale_density,
)
    nspins = rand(Poisson(density.value * scale_density))
    fill(zero(SVector{1, Float64}), nspins)
end

size_scale(::InfiniteWall) = Inf

_geometry_mesh(::InfiniteWall; kwargs...) = [0.]
