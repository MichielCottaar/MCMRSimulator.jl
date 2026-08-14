module GridDispatch

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ...RayGridIntersection: ray_grid_intersections
import ..PhysicalGeometries: Intersection

"""Dispatch a ray through a precomputed grid and stop at the first hit voxel."""
function detect_intersection_grid(
    detect_child,
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N},
    inv_resolution::SVector{N, Float64},
    indices::Array{Vector{Int}, N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{N},
) where {N}
    iszero(destination - start) && return Intersection{N}()
    InternalBoundingBoxes.could_intersect(bounding_box, start, destination) || return Intersection{N}()
    InternalBoundingBoxes.does_intersect(bounding_box, start, destination) || return Intersection{N}()

    lower_bound = InternalBoundingBoxes.lower(bounding_box)
    scaled_start = (start - lower_bound) .* inv_resolution
    scaled_destination = (destination - lower_bound) .* inv_resolution
    found = Intersection{N}()
    entered_grid = false
    grid_size = SVector{N, Int}(size(indices))

    for (voxel, _, _, exit_time, _) in ray_grid_intersections(scaled_start, scaled_destination)
        if any(voxel .< 0) || any(voxel .>= grid_size)
            entered_grid && break
            continue
        end
        entered_grid = true

        for child_index in indices[voxel .+ 1...]
            intersection = detect_child(child_index, previous_hit)
            if intersection.distance < found.distance
                found = intersection
            end
        end
        found.distance <= exit_time && return found
    end
    found
end

end
