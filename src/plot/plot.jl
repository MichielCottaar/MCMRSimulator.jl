using Makie

"""
    color(orient::SpinOrientation; saturation=1.)

Returns a color representing the spin orientation in the transverse (x-y) plane.
Brighter colors have a larger transverse component, so that spins with no transverse component are black.
The actual color encodes the spin orientation.
"""
color(orient::Union{Spin{1}, SpinOrientation}; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))

include("plot_plane.jl")
include("sequence.jl")
include("geometry.jl")
include("off_resonance.jl")
include("snapshot.jl")
include("trajectory.jl")
include("movie.jl")