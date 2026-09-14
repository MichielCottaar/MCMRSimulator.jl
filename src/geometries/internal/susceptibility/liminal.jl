"""Susceptibility state for liminal geometries."""
module Liminal

import StaticArrays: SVector
import ..Grid: susceptibility_off_resonance, off_resonance_gradient

struct LiminalSusceptibility{C}
    children::C
end

function susceptibility_off_resonance(
    susceptibility::LiminalSusceptibility,
    position::SVector{3, Float64},
    previous_hit=nothing,
    isinside=nothing,
)
    if isnothing(previous_hit)
        isnothing(isinside) && return 0.0
        isempty(isinside.inside_of) && return 0.0
        cell_index, offset = first(isinside.inside_of)[1]
        child_previous_hit = nothing
    else
        cell_index, offset = previous_hit.indices[1]
        child_previous_hit = previous_hit
    end
    child = susceptibility.children[cell_index]
    susceptibility_off_resonance(child, position - offset, child_previous_hit)
end

function off_resonance_gradient(susceptibility::LiminalSusceptibility, B0)
    maximum(
        off_resonance_gradient(child, B0)
        for child in susceptibility.children;
        init=0.0,
    )
end

end
