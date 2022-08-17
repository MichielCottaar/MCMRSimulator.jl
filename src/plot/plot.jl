import Makie: Makie, @lift

"""
    color(orient::SpinOrientation; saturation=1.)

Returns a color representing the spin orientation in the transverse (x-y) plane.
Brighter colors have a larger transverse component, so that spins with no transverse component are black.
The actual color encodes the spin orientation.
"""
color(orient::Union{Spin, SpinOrientation}; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))

include("plot_plane.jl")
include("sequence.jl")
include("snapshot.jl")