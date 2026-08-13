module Indices
import StaticArrays: SVector

"""
The index of an base obstruction with the `FixedGeometry`.

This is used internally to identify the obstruction that a particle bumped into.
"""
struct ObstructionIndex{N}
    indices::SVector{N, Int}
end

ObstructionIndex() = ObstructionIndex{0}(SVector{0, Int}())

"""
    get_obstruction(geometry, obstruction_index)

Gets the obstruction identified from `obstruction_index` from within the `PhysicalGeometry`.
"""
function get_obstruction() end

end
