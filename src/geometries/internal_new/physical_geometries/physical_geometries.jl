module PhysicalGeometries

"""
    PhysicalGeometry{N}

Represents a physical geometry in N-dimensional space. This is an abstract type that can be extended to define specific geometries and their properties.
"""
abstract type PhysicalGeometry{N} end

include("groups.jl")
include("transformations.jl")

end
