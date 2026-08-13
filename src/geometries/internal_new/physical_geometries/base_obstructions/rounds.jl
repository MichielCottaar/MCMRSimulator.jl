struct Round{N} <: BaseObstruction{N}
    radius::Float64
end

const InfiniteCylinder = Round{2}
const Sphere = Round{3}

InternalBoundingBox(round::Round{N}) where {N} = InternalBoundingBox{N}(round.radius)

isinside(round::Round{N}, position::SVector{N, Float64}) where {N} =
    sum(position .* position) < round.radius^2

function detect_intersection(
    round::Round{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    obstruction_index::ObstructionIndex=ObstructionIndex(),
) where {N}
    inside = isinside(round, start)
    difference = destination - start
    a = sum(difference .* difference)
    b = sum(2 .* start .* difference)
    c = sum(start .* start)
    determinant = b * b - 4 * a * (c - round.radius^2)
    determinant < 0 && return Intersection{N}()

    solution = (inside ? -b + sqrt(determinant) : -b - sqrt(determinant)) / (2 * a)
    (solution <= 0 || solution > 1) && return Intersection{N}()

    normal = (solution .* destination .+ (1 - solution) .* start) ./ round.radius
    return Intersection(solution, inside ? -normal : normal, inside, obstruction_index, false)
end
