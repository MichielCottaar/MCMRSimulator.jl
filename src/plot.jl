"""
Defines the plotting functions.

This modules does not actually implement these functions.
Instead they are implemented in the `MakieMCMRSimulator` extension if [`Makie`](https://makie.org) is installed and used.

The exception to this is the [`PlotPlane`](@ref), which is defined here.
"""
module Plot
import CoordinateTransformations
import StaticArrays: SVector, MVector
import ..Spins: Snapshot, orientation, SpinOrientation, position
import ..Methods: get_rotation
import ..Geometries.Internal: ray_grid_intersections

"""
    plot_geometry([plot_plane], geometry)

Plots the geometry in a new plot. 
If `plot_plane` is provided the projection on it will be plotted. Otherwise, the 3-dimensional geometry will be plotted.

This function will only work if [`Makie`](https://makie.org) is installed and used.
"""
function plot_geometry end

"""
    plot_geometry!(axis, [plot_plane], geometry)

Plots the geometry on an existing plot axis.
If `plot_plane` is provided the projection on it will be plotted. Otherwise, the 3-dimensional geometry will be plotted.

This function will only work if [`Makie`](https://makie.org) is installed and used.
"""
function plot_geometry! end

"""
    plot_sequence(sequence)

Plots the sequence in a new plot. 

This function will only work if [`Makie`](https://makie.org) is installed and used.
"""
function plot_sequence end

"""
    plot_sequence!(axis, sequence)

Plots the sequence on an existing plot axis.

This function will only work if [`Makie`](https://makie.org) is installed and used.
"""
function plot_sequence! end

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