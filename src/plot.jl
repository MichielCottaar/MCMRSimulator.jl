"""
Defines the plotting functions.

This modules does not actually implement these functions.
Instead they are implemented in the `MakieMCMRSimulator` extension if [`Makie`](https://makie.org) is installed.
"""
module Plot

"""
    plot_geometry([plot_plane], geometry)

Plots the geometry in a new plot. 
If `plot_plane` is provided the projection on it will be plotted. Otherwise, the 3-dimensional geometry will be plotted.

This function will only work if [`Makie`](https://makie.org) is installed.
"""
function plot_geometry end

"""
    plot_geometry!(axis, [plot_plane], geometry)

Plots the geometry on an existing plot axis.
If `plot_plane` is provided the projection on it will be plotted. Otherwise, the 3-dimensional geometry will be plotted.

This function will only work if [`Makie`](https://makie.org) is installed.
"""
function plot_geometry! end

end