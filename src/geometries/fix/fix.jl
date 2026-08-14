module Fix

include("fix_base_geometry.jl")
include("fix_transformations.jl")
include("fix_properties.jl")

import .FixBaseGeometry: fix_base_geometry
import .FixTransformations: fix_transformations
import .FixProperties: fix_properties

end
