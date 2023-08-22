module Susceptibility
include("base.jl")
include("parent.jl")
include("cylinder.jl")
include("annulus.jl")

import .Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient
import .Parent: FixedSusceptibility, ParentSusceptibility, susceptibility_off_resonance, off_resonance_gradient
import .Cylinder: CylinderSusceptibility
import .Annulus: AnnulusSusceptibility

end