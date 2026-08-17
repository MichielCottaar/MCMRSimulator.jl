struct InfiniteWall <: BaseObstruction{1} end

has_inside(::Type{InfiniteWall}) = false
has_single_inside(::Type{InfiniteWall}) = false

const _negative_wall_normal = SVector(-1.0)
const _positive_wall_normal = SVector(1.0)

InternalBoundingBox(::InfiniteWall) = InternalBoundingBox{1}(SVector(0.0), SVector(0.0))

function detect_intersection(
    ::InfiniteWall,
    start::SVector{1, Float64},
    destination::SVector{1, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
)
    !Base.isempty(previous_hit) && return Intersection{1}()
    origin = start[1]
    destination_position = destination[1]
    origin * destination_position > 0 && return Intersection{1}()
    inside = origin > 0
    return Intersection(
        abs(origin) / abs(origin - destination_position),
        inside ? _positive_wall_normal : _negative_wall_normal,
        inside,
        ObstructionIndex(),
        false,
    )
end

function surface_sampling(
    ::InfiniteWall, density::GeometryLeafProperties,
    bounding_box::InternalBoundingBox{1}, scale_density,
)
    nspins = rand(Poisson(density.value * scale_density))
    fill(zero(SVector{1, Float64}), nspins), fill(SVector{1, Float64}(1.0), nspins)
end

size_scale(::InfiniteWall) = Inf

_geometry_mesh(::InfiniteWall; kwargs...) = [0.]
