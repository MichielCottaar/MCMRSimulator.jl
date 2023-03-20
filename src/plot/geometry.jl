"""
    plot(plot_plane, geometry)
    plot!(plot_plane, geometry)
    plot_geometry(plot_plane, geometry)
    plot_geometry!(plot_plane, geometry)

Plots the intersections of `geometry` in the [`PlotPlane`](@ref).
"""
@Makie.recipe(Plot_Geometry, plot_plane, geometry) do scene
    Makie.Theme(
    )
end

function Makie.plot!(pg::Plot_Geometry)
    plot_plane = pg[1]
    raw_geometry = pg[2]
    geometry = @lift Geometry($raw_geometry)

    to_plot = @lift project_geometry($plot_plane, $geometry)

    on(to_plot) do to_iter
        for (func, args, kwargs) in to_iter
            func(pg, args...; kwargs...)
        end
    end
    to_plot[] = to_plot[]
    pg
end

Makie.plottype(::PlotPlane, ::Obstruction) = Plot_Geometry
Makie.plottype(::PlotPlane, ::AbstractVector{<:Obstruction}) = Plot_Geometry
Makie.plottype(::PlotPlane, ::TransformObstruction) = Plot_Geometry
Makie.plottype(::PlotPlane, ::AbstractVector{<:TransformObstruction}) = Plot_Geometry
Makie.plottype(::PlotPlane, ::Geometry) = Plot_Geometry


function project_geometry(plot_plane::PlotPlane, geometry::Geometry)
    projections = []
    for t in geometry.obstructions
        append!(projections, project_geometry(plot_plane, t))
    end
    projections
end

"""
    project_geometry(plot_plane, transform)

Projects the plane on the intrinsic plane of the obstructions deformed by `transform`.
"""
function project_geometry(plot_plane::PlotPlane, transform::TransformObstruction{N}) where {N}
    center_obstruction_space = plot_plane.transformation(zero(PosVector))

    obstruction_coordinates_in_plot_plane = SVector{N}(map(p->PosVector(plot_plane.transformation(p)) - center_obstruction_space, eachcol(transform.rotation)))

    projections = []
    for (single, obstruction) in zip(transform.transforms, transform.obstructions)
        obstruction_center_in_plot_plane = plot_plane.transformation(transform.rotation * single.shift)
        append!(projections, project_obstruction(obstruction, obstruction_center_in_plot_plane, obstruction_coordinates_in_plot_plane, transform.repeats, (plot_plane.sizex, plot_plane.sizey)))
    end
    projections
end


function project_obstruction(wall::Wall, center::PosVector, obstruction_coordinates_in_plot_plane::SVector{1, PosVector}, repeats::SVector{1, Float}, sizes::Tuple{<:Real, <:Real})
    constant_center = obstruction_coordinates_in_plot_plane[1] ⋅ center
    normal = obstruction_coordinates_in_plot_plane[1][1:2]
    halfs = (sizes[1]/2, sizes[2]/2)

    function get_line(constant::Float)
        if abs(normal[1]) > abs(normal[2])
            return [
                SVector{2, Float}((constant + normal[2] * sizes[2]) / normal[1], -sizes[2]),
                SVector{2, Float}((constant - normal[2] * sizes[2]) / normal[1], sizes[2]),
                SVector{2, Float}(NaN, NaN),
            ]
        else
            return [
                SVector{2, Float}(sizes[1], (constant - normal[1] * sizes[1]) / normal[2]),
                SVector{2, Float}(-sizes[1], (constant + normal[1] * sizes[1]) / normal[2]),
                SVector{2, Float}(NaN, NaN),
            ]
        end
    end
    if isfinite(repeats[1])
        repeat_constant_shift = repeats[1] * norm(normal)
        closest_center = mod(constant_center, repeat_constant_shift)
        max_size = abs.(normal) ⋅ halfs
        nshift = div(max_size + abs(closest_center), repeat_constant_shift, RoundUp)
    
        line = SVector{2, Float}[]
        for shift in -nshift:nshift
            append!(line, get_line(shift * repeat_constant_shift + closest_center))
        end
    else
        line = get_line(constant_center)[1:2]
    end
    [(Makie.lines!, (cut_line(line, halfs), ), Dict())]
end

function project_obstruction(obstruction::Union{Cylinder, Spiral}, center_vec::PosVector, obstruction_coordinates_in_plot_plane::SVector{2, PosVector}, repeats::SVector{2, Float}, sizes::Tuple{<:Real, <:Real})
    (dirx, diry) = obstruction_coordinates_in_plot_plane
    normal = cross(dirx, diry)
    if normal[3] == 0
        error("Can not plot cylinders or spirals perfectly aligned with plotting plane")
    end
    halfs = (sizes[1]/2, sizes[2]/2)
    center = (center_vec .- normal .* (center_vec[3] / normal[3]))[1:2]
    dirx = (dirx .- normal .* (dirx[3] / normal[3]))[1:2]
    diry = (diry .- normal .* (diry[3] / normal[3]))[1:2]

    projection_size = normal[3] * normal[3]

    theta_normal = atan(normal[2], normal[1])
    function relative_radius(theta)
        ct_rot = cos(theta - theta_normal)
        sqrt(1 / (1 - ct_rot * ct_rot * (1 - projection_size)))
    end
    if isa(obstruction, Cylinder)
        theta = (0:0.02:2) * π
        radius = obstruction.radius .* relative_radius.(theta)
    elseif isa(obstruction, Spiral)
        theta = obstruction.theta0:0.02:obstruction.theta_end
        radius = @. (theta - obstruction.theta0) * obstruction.thickness / 2π + obstruction.inner
    end

    points = map((r, t) -> SVector{2, Float}(r .* cos(t), r * sin(t)), radius, theta)
    push!(points, SVector{2, Float}(NaN, NaN))

    function get_line(center :: SVector{2, Float})
        map(p->SVector{2, Float}(p[1] + center[1], p[2] + center[2]), points)
    end

    if all(isfinite.(repeats))
        stepx = dirx * repeats[1]
        stepy = diry * repeats[2]

        toshift = round.([stepx stepy] \ center)
        closest_center = center - (toshift[1] * stepx + toshift[2] * stepy)

        outer_radius = isa(obstruction, Cylinder) ? obstruction.radius : obstruction.outer
        nshiftx = div(abs.(dirx) ⋅ (halfs .+ abs.(closest_center)) + projection_size * outer_radius, repeats[1], RoundUp)
        nshifty = div(abs.(diry) ⋅ (halfs .+ abs.(closest_center)) + projection_size * outer_radius, repeats[2], RoundUp)
        line = SVector{2, Float}[]
        for i in -nshiftx:nshiftx
            int_pos = closest_center .+ i .* stepx
            for j in -nshifty:nshifty
                pos = SVector{2}(int_pos .+ j .* stepy)
                append!(line, get_line(pos))
            end
        end
    else
        line = get_line(SVector{2}(center))
    end
    res = [(Makie.lines!, (cut_line(line, halfs), ), Dict())]
    if isa(obstruction, Spiral)
        if obstruction.inner_cylinder
            append!(res, project_obstruction(obstruction.equivalent_annulus.inner, center_vec, obstruction_coordinates_in_plot_plane, repeats, sizes))
        end
        if obstruction.outer_cylinder
            append!(res, project_obstruction(obstruction.equivalent_annulus.outer, center_vec, obstruction_coordinates_in_plot_plane, repeats, sizes))
        end
    end
    res
end

function project_obstruction(annulus::Annulus, args...)
    [
        project_obstruction(annulus.inner, args...)[1],
        project_obstruction(annulus.outer, args...)[1],
    ]
end

function cut_line(old_line::AbstractVector{SVector{2, Float}}, sizes)
    _inside(point::SVector{2, Float}) = abs(point[1]) <= sizes[1] && abs(point[2]) <= sizes[2]
    new_line = SVector{2, Float}[]
    prev_point = SVector{2, Float}(NaN, NaN)

    for new_point in old_line
        if !all(isfinite.(new_point))
            if length(new_line) > 0 && all(isfinite.(new_line[end]))
                push!(new_line, new_point)
            end
        else
            if !all(isfinite.(prev_point))
                if _inside(new_point)
                    push!(new_line, new_point)
                end
            else
                # both prev_point and new_point are finite
                add_overlap!(new_line, prev_point, new_point, sizes)
            end
        end
        prev_point = new_point
    end
    map(p->Makie.Point2(p[1], p[2]), new_line)
end

function add_overlap!(new_line, prev_point, new_point, sizes)
    _inside(point::SVector{2, Float}) = abs(point[1]) <= sizes[1] && abs(point[2]) <= sizes[2]

    next_inside = _inside(new_point)
    ndiscover = 2 - _inside(prev_point) - next_inside
    if ndiscover == 0
        push!(new_line, new_point)
        return
    end

    # line is at least partially outside
    direction = new_point - prev_point
    dist_new_sq = sum(d->d*d, direction)
    n = norm(direction)
    normal = SVector{2, Float}(-direction[2]/n, direction[1]/n)
    constant = normal ⋅ prev_point

    discovered = SVector{2, Float}[]
    function add_valid!(point::SVector{2, Float})
        distsq = (point - prev_point) ⋅ direction
        if distsq >= 0 && distsq <= dist_new_sq
            push!(discovered, point)
        end
    end
    x1 = (constant - sizes[2] * normal[2]) / normal[1]
    if abs(x1) <= sizes[1]
        add_valid!(SVector{2, Float}(x1, sizes[2]))
    end
    x2 = (constant + sizes[2] * normal[2]) / normal[1]
    if abs(x2) <= sizes[1]
        add_valid!(SVector{2, Float}(x2, -sizes[2]))
    end

    y1 = (constant - sizes[1] * normal[1]) / normal[2]
    if abs(y1) < sizes[2]
        add_valid!(SVector{2, Float}(sizes[1], y1))
    end
    y2 = (constant + sizes[1] * normal[1]) / normal[2]
    if abs(y2) < sizes[2]
        add_valid!(SVector{2, Float}(-sizes[1], y2))
    end

    if length(discovered) == 0 && ndiscover == 2
        # line is fully outside of the box
        return
    end
    if length(discovered) != ndiscover
        # bad point at edge; let's just ignore it
        return
    end

    if ndiscover == 2
        if length(new_line) > 0 && all(isfinite.(new_line[end]))
            push!(new_line, SVector{2, Float}(NaN, NaN))
        end
        if (discovered[1] - prev_point) ⋅ direction > (discovered[2] - prev_point) ⋅ direction
            push!(new_line, discovered[2])
            push!(new_line, discovered[1])
        else
            push!(new_line, discovered[1])
            push!(new_line, discovered[2])
        end
    else
        @assert ndiscover == 1
        if next_inside && length(new_line) > 0 && all(isfinite.(new_line[end]))
            push!(new_line, SVector{2, Float}(NaN, NaN))
        end
        push!(new_line, discovered[1])
        if next_inside
            push!(new_line, new_point)
        end
    end
end