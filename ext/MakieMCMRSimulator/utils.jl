module Utils
import Makie
import Colors
import MCMRSimulator.Spins: Spin, SpinOrientation, phase, transverse
import MCMRSimulator: FixedGeometry, Snapshot, ObstructionGroup
import MCMRSimulator.Plot: PlotPlane


"""
    color(orient::SpinOrientation; saturation=1.)
    color(spin::Spin; saturation=1., sequence=1)

Returns a color representing the spin orientation in the transverse (x-y) plane.
Brighter colors have a larger transverse component, so that spins with no transverse component are black.
The actual color encodes the spin orientation.
"""
color(orient::SpinOrientation; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))
color(spin::Spin{N}; sequence=1, kwargs...) where {N} = color(spin.orientations[sequence]; kwargs...)


const Projectable = Union{FixedGeometry, ObstructionGroup, Vector{<:ObstructionGroup}, Snapshot, Vector{<:Snapshot}}

Makie.args_preferred_axis(::Projectable) = Makie.LScene
Makie.args_preferred_axis(::PlotPlane) = Makie.Axis

end