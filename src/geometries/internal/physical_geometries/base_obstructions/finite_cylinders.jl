import LinearAlgebra: cross, norm, ⋅

struct FiniteCylinder <: BaseObstruction{3}
    first::SVector{3, Float64}
    axis::SVector{3, Float64}
    length::Float64
    radius_first::Float64
    radius_second::Float64
    caps_are_gaps::Bool
end

function FiniteCylinder(
    first::AbstractVector, second::AbstractVector,
    radius::Real;
    radius_second::Real=radius,
    caps_are_gaps=false,
)
    first = SVector{3, Float64}(first)
    second = SVector{3, Float64}(second)
    displacement = second - first
    length = norm(displacement)
    length > 0 || throw(ArgumentError("finite-cylinder endpoints must be distinct"))
    radius > 0 || throw(ArgumentError("finite-cylinder radius must be positive"))
    radius_second > 0 || throw(ArgumentError("finite-cylinder radius must be positive"))
    caps_are_gaps isa Bool || throw(ArgumentError("caps_are_gaps must be a Bool"))
    FiniteCylinder(first, displacement / length, length, Float64(radius), Float64(radius_second), caps_are_gaps)
end

FiniteCylinder(
    first::AbstractVector, second::AbstractVector,
    radius_first::Real, radius_second::Real;
    caps_are_gaps=false,
) = FiniteCylinder(first, second, radius_first; radius_second, caps_are_gaps)

FiniteCylinder(
    first::AbstractVector, second::AbstractVector;
    radius,
    radius_second=radius,
    caps_are_gaps=false,
) = FiniteCylinder(first, second, radius, radius_second; caps_are_gaps)

has_inside(::Type{FiniteCylinder}) = true
has_single_inside(::Type{FiniteCylinder}) = true
intersection_type(::Type{FiniteCylinder}) = Tuple{Int}

_finite_cylinder_second(cylinder::FiniteCylinder) =
    cylinder.first + cylinder.length * cylinder.axis

function _finite_cylinder_radius(cylinder::FiniteCylinder, axial_position)
    fraction = axial_position / cylinder.length
    (1 - fraction) * cylinder.radius_first + fraction * cylinder.radius_second
end

function isinside_single(
    cylinder::FiniteCylinder,
    position::SVector{3, Float64},
    previous_intersection=nothing,
)
    !isnothing(previous_intersection) && return previous_intersection[2]
    relative = position - cylinder.first
    axial = relative ⋅ cylinder.axis
    (0 < axial < cylinder.length) || return false
    radial = relative - axial * cylinder.axis
    radial ⋅ radial < _finite_cylinder_radius(cylinder, axial)^2
end

function InternalBoundingBox(cylinder::FiniteCylinder)
    center = cylinder.first + cylinder.length / 2 * cylinder.axis
    half_size = cylinder.length / 2 .* abs.(cylinder.axis) .+
        max(cylinder.radius_first, cylinder.radius_second) .* sqrt.(1 .- cylinder.axis .* cylinder.axis)
    InternalBoundingBox(half_size, center)
end

size_scale(cylinder::FiniteCylinder) = min(cylinder.length, min(cylinder.radius_first, cylinder.radius_second))

function _finite_cylinder_side_distance(cylinder, axial, radial_distance)
    radius_difference = cylinder.radius_second - cylinder.radius_first
    radius_slope = radius_difference / cylinder.length
    direction = SVector(1.0, radius_slope)
    point = SVector(axial, radial_distance - cylinder.radius_first)
    fraction = clamp((point ⋅ direction) / (direction ⋅ direction), 0., 1.)
    norm(point - fraction .* direction)
end

function distance_to_surface(cylinder::FiniteCylinder, position::SVector{3, Float64})
    relative = position - cylinder.first
    axial = relative ⋅ cylinder.axis
    radial = relative - axial * cylinder.axis
    radial_distance = norm(radial)
    side_distance = _finite_cylinder_side_distance(cylinder, axial, radial_distance)
    first_cap_distance = axial <= 0 ?
        hypot(-axial, max(radial_distance - cylinder.radius_first, 0.)) :
        hypot(axial, max(radial_distance - cylinder.radius_first, 0.))
    second_axial = axial - cylinder.length
    second_cap_distance = hypot(second_axial, max(radial_distance - cylinder.radius_second, 0.))
    min(side_distance, first_cap_distance, second_cap_distance)
end

function _finite_cylinder_candidate(current, index, distance)
    (0 < distance <= 1 && (isnothing(current) || distance < current[2])) ?
        (index, distance) : current
end

function _finite_cylinder_side_candidate(current, cylinder, axial_start, axial_displacement, distance)
    axial = axial_start + distance * axial_displacement
    0 <= axial <= cylinder.length || return current
    _finite_cylinder_candidate(current, 1, distance)
end

function find_intersection(
    cylinder::FiniteCylinder,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit=nothing,
)
    previous = !isnothing(previous_hit)
    inside = previous ? previous_hit[2] : isinside_single(cylinder, start)
    !inside && previous && return nothing

    displacement = destination - start
    relative = start - cylinder.first
    axial_start = relative ⋅ cylinder.axis
    axial_displacement = displacement ⋅ cylinder.axis
    radial_start = relative - axial_start * cylinder.axis
    radial_displacement = displacement - axial_displacement * cylinder.axis
    radius_start = _finite_cylinder_radius(cylinder, axial_start)
    radius_slope = (cylinder.radius_second - cylinder.radius_first) / cylinder.length * axial_displacement

    a = radial_displacement ⋅ radial_displacement - radius_slope^2
    b = 2 * (radial_start ⋅ radial_displacement - radius_start * radius_slope)
    c = radial_start ⋅ radial_start - radius_start^2
    best = nothing
    if abs(a) < 1e-12
        abs(b) > 1e-12 && (best = _finite_cylinder_side_candidate(best, cylinder, axial_start, axial_displacement, -c / b))
    else
        determinant = b^2 - 4 * a * c
        if determinant >= 0
            root = sqrt(determinant)
            best = _finite_cylinder_side_candidate(best, cylinder, axial_start, axial_displacement, (-b - root) / (2 * a))
            best = _finite_cylinder_side_candidate(best, cylinder, axial_start, axial_displacement, (-b + root) / (2 * a))
        end
    end

    if abs(axial_displacement) > 1e-12
        first_cap = -axial_start / axial_displacement
        second_cap = (cylinder.length - axial_start) / axial_displacement
        radial_first = radial_start + first_cap * radial_displacement
        radial_second = radial_start + second_cap * radial_displacement
        if radial_first ⋅ radial_first <= cylinder.radius_first^2
            best = _finite_cylinder_candidate(best, 2, first_cap)
        end
        if radial_second ⋅ radial_second <= cylinder.radius_second^2
            best = _finite_cylinder_candidate(best, 3, second_cap)
        end
    end
    isnothing(best) ? nothing : (best[1], inside, best[2])
end

function get_intersection_params(
    cylinder::FiniteCylinder,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    intersection::Tuple,
    isinside=nothing,
)
    index, inside, distance = intersection
    position = (1 - distance) .* start + distance .* destination
    outward_normal = if index == 1
        relative = position - cylinder.first
        axial = relative ⋅ cylinder.axis
        radial = relative - axial * cylinder.axis
        radial_normal = radial / norm(radial)
        slope = (cylinder.radius_second - cylinder.radius_first) / cylinder.length
        (radial_normal - slope * cylinder.axis) / sqrt(1 + slope^2)
    elseif index == 2
        -cylinder.axis
    elseif index == 3
        cylinder.axis
    else
        throw(ArgumentError("invalid finite-cylinder intersection index"))
    end
    IntersectionParams{3}(inside, inside ? -outward_normal : outward_normal, index != 1 && cylinder.caps_are_gaps)
end

function _finite_cylinder_basis(axis)
    reference = abs(axis[1]) < 0.9 ? SVector(1., 0., 0.) : SVector(0., 1., 0.)
    first = cross(axis, reference)
    first / norm(first), cross(axis, first / norm(first))
end

function surface_sampling(
    cylinder::FiniteCylinder, density::GeometryLeafProperties,
    scale_density,
)
    radius_difference = cylinder.radius_second - cylinder.radius_first
    slant = hypot(cylinder.length, radius_difference)
    side_area = π * (cylinder.radius_first + cylinder.radius_second) * slant
    cap_area = cylinder.caps_are_gaps ? 0. : π * (cylinder.radius_first^2 + cylinder.radius_second^2)
    total_area = side_area + cap_area
    nspins = rand(Poisson(total_area * density.value * scale_density))
    first_basis, second_basis = _finite_cylinder_basis(cylinder.axis)
    positions = Vector{SVector{3, Float64}}()
    for _ in 1:nspins
        if rand() * total_area < side_area
            axial = rand() * cylinder.length
            radius = _finite_cylinder_radius(cylinder, axial)
            theta = rand() * 2π
            push!(positions, cylinder.first + axial * cylinder.axis + radius * (cos(theta) * first_basis + sin(theta) * second_basis))
        else
            at_second = rand(Bool)
            radius = at_second ? cylinder.radius_second : cylinder.radius_first
            theta = rand() * 2π
            radial = sqrt(rand()) * radius * (cos(theta) * first_basis + sin(theta) * second_basis)
            push!(positions, cylinder.first + (at_second ? cylinder.length : 0.) * cylinder.axis + radial)
        end
    end
    positions
end

function _geometry_mesh(cylinder::FiniteCylinder; nsamples=100, kwargs...)
    nsamples >= 3 || throw(ArgumentError("nsamples must be at least 3"))
    first_basis, second_basis = _finite_cylinder_basis(cylinder.axis)
    first_circle = [cylinder.first + cylinder.radius_first * (cos(2π * i / nsamples) * first_basis + sin(2π * i / nsamples) * second_basis) for i in 0:(nsamples - 1)]
    second_circle = [cylinder.first + cylinder.length * cylinder.axis + cylinder.radius_second * (cos(2π * i / nsamples) * first_basis + sin(2π * i / nsamples) * second_basis) for i in 0:(nsamples - 1)]
    vertices = vcat(first_circle, second_circle, [cylinder.first, _finite_cylinder_second(cylinder)])
    triangles = SVector{3, Int}[]
    for i in 1:nsamples
        next = mod1(i + 1, nsamples)
        push!(triangles, SVector(i, next, nsamples + i))
        push!(triangles, SVector(next, nsamples + next, nsamples + i))
        push!(triangles, SVector(2 * nsamples + 1, next, i))
        push!(triangles, SVector(2 * nsamples + 2, nsamples + i, nsamples + next))
    end
    [_mesh_result(vertices, triangles)]
end
