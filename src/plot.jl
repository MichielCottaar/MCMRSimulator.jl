"""
This module contains the base plotting functions.

If [`Makie`](https://makie.org) is not installed, none of the plotting functions will work.
If a Makie backend is installed and importend, concrete methods will be added to these functions by the MakieMCMRSimulator extension.
This makes Makie an optional dependency of MCMRSimulator, which will only be required if you want to use the MCMRSimulator plotting capabilities.

In addition to these empty plotting functions, this module defines the [`PlotPlane`](@ref) for any 2D-projections and helper functions to project onto this plane.
"""
module Plot
import CoordinateTransformations
import StaticArrays: SVector, MVector
import ...Spins: Snapshot, orientation, SpinOrientation, position
import ...Methods: get_rotation
import ...Geometries.Internal: ray_grid_intersections


function plot_geometry end
function plot_geometry! end
"""
    plot([plot_plane,] geometry; kwargs...)
    plot!([scene,] [plot_plane,] geometry; kwargs...)
    plot_geometry([plot_plane,] geometry; kwargs...)
    plot_geometry!([scene,] [plot_plane,] geometry; kwargs...)

Plots a given geometry. 
If a [`PlotPlane`](@ref) is provided the 2D projection of the geometry onto this plane is plotted.
Otherwise, the geometry is plotted in 3D.

If you want to overlay the off-resonance field, call [`plot_off_resonance`](@ref) first before calling this function.

This function will only work if a [`Makie`](https://makie.org) backend is imported.
"""
plot_geometry, plot_geometry!

function plot_off_resonance end
function plot_off_resonance! end
"""
    plot_off_resonance([plot_plane,] geometry; kwargs...)
    plot_off_resonance!([scene,] [plot_plane,] geometry; kwargs...)

Plots the off-resonance field generated by a given geometry. 
If a [`PlotPlane`](@ref) is provided the 2D projection of the geometry onto this plane is plotted.
Otherwise, the geometry is plotted in 3D.

This function will only work if a [`Makie`](https://makie.org) backend is imported.
"""
plot_off_resonance, plot_off_resonance!


function plot_snapshot end
function plot_snapshot! end
"""
    plot([plot_plane], snapshot; kind=:scatter, kwargs...)
    plot!([scene,] [plot_plane], snapshot; kind=:scatter, kwargs...)
    plot_snapshot([plot_plane], snapshot; kind=:scatter, kwargs...)
    plot_snapshot!([scene,] [plot_plane], snapshot; kind=:scatter, kwargs...)

Plots a [`Snapshot`](@ref) in a new plot. 

The spins are plotted in 2D projected onto the [`PlotPlane`](@ref) if one is provided.
Otherwise, the spins are plotted in 3D (does not work for `kind=:image`).

There are three kinds of snapshot plots available:
## Scatter plot
Default (or set using `kind=:scatter`). Each spin is plotted as a point with the colour set by the transverse magnetisation.
Additional keywords are passed on to `Makie.scatter`.

## Dyad plot
Set using `kind=:dyad`. Each spin is plotted as a dyad. For a 2D dyad the orienation is set by the transverse magnetisation.
For a 3D dyad the full magnetisation is used to set the orienation.
The length of the dyads can be controlled using `dyad_length` (0.1 by default).
Additional keywords are passed on to `Makie.arrows`.

## Image
Set using `kind=:image`. The average magnetisation is plotted across the [`PlotPlane`](@ref). 
The colour in each pixel is set by the average transverse magnetisation of the local spins.
Additional keywords are passed on to `Makie.image!`.
This plot will not work in 3D (i.e., a [`PlotPlane`](@ref) is required).

This function will only work if [`Makie`](https://makie.org) is installed and imported.
"""
plot_snapshot, plot_snapshot!

function plot_trajectory end
function plot_trajectory! end
"""
    plot([plot_plane], snapshots; kwargs...)
    plot!([scene,] [plot_plane], snapshots; kwargs...)
    plot_trajectory([plot_plane], snapshots; kwargs...)
    plot_trajectory!([scene,] [plot_plane], snapshots; kwargs...)

Plots the spin trajectory in a vector of [`Snapshot`](@ref) on an existing plot scene.

The spins are plotted in 2D projected onto the [`PlotPlane`](@ref) if one is provided.
Otherwise, the spins are plotted in 3D.
At each location along the trajectory, the colour is set by the transverse magnetisation.
Additional keywords are passed on to `Makie.lines!`.

This function will only work if [`Makie`](https://makie.org) is installed and imported.
"""
plot_trajectory, plot_trajectory!

function plot_sequence end
function plot_sequence! end
"""
    plot(sequence; kwargs...)
    plot!([scene,] sequence; kwargs...)
    plot_sequence(sequence; kwargs...)
    plot_sequence!([scene,] sequence; kwargs...)

Plot the sequence diagram.

This function will only work if [`Makie`](https://makie.org) is installed and imported.
"""
plot_sequence, plot_sequence!

"""
    simulator_movie(filename, simulation, times, size; resolution=(1600, 800), trajectory_init=30, signal_init=10000, framerate=50, plane_orientation=:z, kwargs...)

Writes a movie of the [`Simulation`](@ref) to the given `filename`.

Each frame of the movie shows the [`Snapshot`](@ref) at given `times`.
`size` is a tuple with the size of the plotted region in the x- and y-direction.
If there is a repeating geometry, then it is strongly recommended to use the size of the repeat for `size`.

Keyword arguments:
- `resolution`: pixel resolution of each frame in the movie.
- `trajectory_init`: how many spins to plot on each frame.
- `signal_init`: how many spins to use to evaluate the signal evolution.
- `framerate`: wait time between each subsequent frame in the movie.
- `plane_orientation`: orienation of the plane on which the spins are projected (see [`PlotPlane`](@ref)).

Additional keyword arguments are passed on to [`plot_snapshot!`](@ref).

This function will only work if [`Makie`](https://makie.org) is installed and imported.
"""
function simulator_movie end


"""
Defines a finite plane in the 3D space used for plotting.

# Constructor
    PlotPlane(normal=[0, 0, 1], position=[0, 0, 0]; size=10., sizex=<size>, sizey=<size>, ngrid=100)

Arguments:
- `normal`: length-3 vector with the orientation perpendicular to the plane (default to z-direction).
- `position`: position of plane as a length-3 vector (defaults to origin).
- `sizex`: size of the plane in the x-direction (before rotating to `normal`).
- `sizey`: size of the plane in the y-direction (before rotating to `normal`).
- `size`: set default value for `sizex` and `sizey`.
- `ngrid`: number of grid elements to split the plane up into for plotting.
"""
struct PlotPlane
    transformation :: CoordinateTransformations.Transformation
    sizex :: Real
    sizey :: Real
end


function PlotPlane(
    normal=:z, 
    position :: AbstractVector{<:Real}=zero(SVector{3, Float64});
    sizex=nothing, sizey=nothing, size=10.,
)
    if isnothing(sizex)
        sizex = size
    end
    if isnothing(sizey)
        sizey = size
    end
    rotation = get_rotation(normal, 3; reference_dimension=:z)
    position = SVector{3}(position)
    transform = CoordinateTransformations.AffineMap(rotation, position)
    PlotPlane(CoordinateTransformations.inv(transform), sizex, sizey)
end


function project(pp::PlotPlane, pos::SVector{3, Float64})
    base = pp.transformation(pos)
    sizex, sizey = pp.sizex, pp.sizey
    correct = [isfinite(sizex) ? sizex : 0., isfinite(sizey) ? sizey : 0., 0.] / 2.
    mod.(base .+ correct, (sizex, sizey, Inf)) .- correct
end

function project_trajectory(pp::PlotPlane, pos::AbstractVector{<:SVector{3, Float64}})
    transformed = pp.transformation.(pos)
    pos2D = map(p -> SVector{3, Float64}(
        isfinite(pp.sizex) ? p[1] / pp.sizex + 0.5 : 0.5,
        isfinite(pp.sizey) ? p[2] / pp.sizey + 0.5 : 0.5,
        0.5
    ), transformed)
    function to_pos!(time::Float64, position::SVector{3, Float64})
        idx1 = Int(floor(time))
        idx2 = Int(ceil(time))
        if idx1 == idx2
            pos_ref = pos2D[idx1]
        else
            pos_ref = pos2D[idx2] .* (time - idx1) .+ pos2D[idx1] .* (idx2 - time)
        end
        push!(all_pos, SVector{2}(
            isfinite(pp.sizex) ? pp.sizex * (position[1] - 0.5) : pos_ref[1],
            isfinite(pp.sizey) ? pp.sizey * (position[2] - 0.5) : pos_ref[2],
        ))
        push!(times, time)
    end
    empty = SVector{2, Float64}(NaN, NaN)
    all_pos = SVector{2, Float64}[]
    times = Float64[]
    for (pos_index, origin, destination) in zip(1:length(pos2D)-1, pos2D[1:end-1], pos2D[2:end])
        for (_, t1, p1, t2, p2) in  ray_grid_intersections(origin, destination)
            if pos_index == 1 && t1 == zero(t1)
                to_pos!(pos_index + t1, p1)
            end
            if t1 == zero(t1) && t2 == one(t2)
                to_pos!(pos_index + t2, p2)
            elseif t1 == zero(t1)
                to_pos!(pos_index + t2, p2)
                push!(all_pos, empty)
                push!(times, NaN)
            elseif t2 == one(t2)
                to_pos!(pos_index + t1, p1)
                to_pos!(pos_index + t2, p2)
            else
                to_pos!(pos_index + t1, p1)
                to_pos!(pos_index + t2, p2)
                push!(all_pos, empty)
                push!(times, NaN)
            end
        end
    end
    return all_pos, times
end


"""
    project_on_grid(plot_plane, snap, ngrid)

Spins from the [`Snapshot`](@ref) are projected onto the grid defined by [`PlotPlane`](@ref) in two ways:
- along the normal spins are projected onto the plane from infinitely far (TODO: give finite extent)
- in the other directions any spins are projected onto the plane using mod(position[1], `sizex`) and mod(position[2], `sizey`).
    This assumes that the geometry and field repeats itself ad infinitum beyond the `PlotPlane` (TODO: allow this assumption to be turned off).
In effect, this means that all spins are projected onto the `PlotPlane`. The average spin orientation in each grid cell is returned.
"""
function project_on_grid(pp::PlotPlane, snap::Snapshot{1}, ngrid::Int)
    positions = [project(pp, position(spin)) for spin in snap]
    xrange = extrema(v->v[1], positions)
    yrange = extrema(v->v[2], positions)

    res = zeros(MVector{3}, ngrid, ngrid)
    hits = zeros(Int, ngrid, ngrid)
    for (spin, pos) in zip(snap, positions)
        relx = (pos[1] - xrange[1])/(xrange[2] - xrange[1])
        rely = (pos[2] - yrange[1])/(yrange[2] - yrange[1])
        x_index = min.(Int(floor(relx * (ngrid - 1))), ngrid - 2) + 1
        y_index = min.(Int(floor(rely * (ngrid - 1))), ngrid - 2) + 1
        vs = orientation(spin)
        for (xi, yi) in [
            (x_index, y_index),
            (x_index+1, y_index),
            (x_index, y_index+1),
            (x_index+1, y_index+1),
        ]
            res[xi, yi] += vs
            hits[xi, yi] += 1
        end
    end
    (range(xrange..., ngrid), range(yrange..., ngrid), SpinOrientation.(res ./ hits))
end
end