module BaseObstructions

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices
import ...InternalBoundingBoxes: InternalBoundingBox
import ...Indices: ObstructionIndex

abstract type BaseObstruction{N} <: PhysicalGeometry{N} end

include("infinite_walls.jl")
include("rounds.jl")
include("overlapping_rounds.jl")
include("triangles.jl")

end
