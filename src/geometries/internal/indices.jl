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
    add_index(index/intersection, new_index::Int)

Add an integer index to the front of the index/intersection. A new object is returned.
"""
add_index(oi::ObstructionIndex{N}, new_index::Int) where {N} = ObstructionIndex{N+1}(SVector{N+1, Int}([new_index, oi.indices...]))

"""
    remove_index(index/intersection)

Returns the first integer index and a new index/intersection with that index removed.
"""
remove_index(oi::ObstructionIndex{N}) where {N} = (oi.indices[1], ObstructionIndex{N-1}(SVector{N-1, Int}([oi.indices[2:end]...])))
remove_index(::ObstructionIndex{0}) = error("Cannot retreive index from empty `ObstructionIntersection`.")

"""
    get_obstruction(geometry, obstruction_index)

Gets the obstruction identified from `obstruction_index` from within the `PhysicalGeometry`.
"""
function get_obstruction() end

end
