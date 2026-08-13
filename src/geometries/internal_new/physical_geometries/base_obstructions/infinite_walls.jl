struct InfiniteWall <: BaseObstruction{1} end

const _negative_wall_normal = SVector(-1.0)
const _positive_wall_normal = SVector(1.0)

InternalBoundingBox(::InfiniteWall) = InternalBoundingBox{1}(SVector(0.0), SVector(0.0))

function detect_intersection(
    ::InfiniteWall,
    start::SVector{1, Float64},
    destination::SVector{1, Float64},
    obstruction_index::ObstructionIndex=ObstructionIndex(),
)
    origin = start[1]
    destination_position = destination[1]
    origin * destination_position > 0 && return Intersection{1}()
    inside = origin > 0
    return Intersection(
        abs(origin) / abs(origin - destination_position),
        inside ? _positive_wall_normal : _negative_wall_normal,
        inside,
        obstruction_index,
        false,
    )
end
