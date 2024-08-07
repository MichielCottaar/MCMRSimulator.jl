module Susceptibility
include("base.jl")
include("parent.jl")
include("cylinder.jl")
include("annulus.jl")
include("triangle.jl")

import .Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient
import .Parent: FixedSusceptibility, ParentSusceptibility, susceptibility_off_resonance, off_resonance_gradient, ShiftedSusceptibility
import .Cylinder: CylinderSusceptibility
import .Annulus: AnnulusSusceptibility
import .Triangle: TriangleSusceptibility

end