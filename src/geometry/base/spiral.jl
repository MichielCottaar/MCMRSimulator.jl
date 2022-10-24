"""
    Spiral(inner, outer; theta0=0., thickness=0.014)

Creates a spiral starting from `inner` to `outer` radius.
The distance between subsequent wraps is `thickness` micrometer.
The angle of the start of the wraps is `theta0`.
"""
struct Spiral <: BaseObstruction{2}
    inner :: Float
    outer :: Float
    theta0 :: Float
    theta_end :: Float
    thickness :: Float
    closed :: Bool
    inner_cylinder :: Bool
    outer_cylinder :: Bool
    equivalent_annulus :: Annulus
    properties :: ObstructionProperties
end

function Spiral(inner, outer; theta0=zero(Float), thickness=0.014, myelin=false, chi_I=-0.1, chi_A=-0.1, inner_cylinder=true, outer_cylinder=false, kwargs...)
    @assert outer > inner
    closed = (outer - inner) > thickness
    theta_end = theta0 + Float(2π) * (outer - inner) / thickness
    Spiral(Float(inner), Float(outer), Float(theta0), Float(theta_end), Float(thickness), closed, inner_cylinder, outer_cylinder, Annulus(inner, outer, myelin=myelin, chi_I=chi_I, chi_A=chi_A), ObstructionProperties(kwargs...))
end

function isinside(s::Spiral, pos::SVector{2, Float})
    rsq = (pos[1] * pos[1]) + (pos[2] * pos[2])
    if !s.closed
        return 0
    elseif rsq > s.outer * s.outer
        return 0
    elseif rsq < s.inner * s.inner
        return 2
    end

    half_inner = s.inner + s.thickness
    half_outer = s.outer - s.thickness
    if rsq < half_outer * half_outer && rsq > half_inner * half_inner
        return 1
    end
    theta = atan(pos[2], pos[1])
    if rsq > half_outer * half_outer
        dtheta = mod((s.theta_end - theta) / 2π, 1)
        local_r = s.outer - dtheta * s.thickness
        return rsq > (local_r * local_r) ? 0 : 1
    else
        dtheta = mod((theta - s.theta0) / 2π, 1)
        local_r = s.inner + dtheta * s.thickness
        return rsq > (local_r * local_r) ? 1 : 2
    end
end


"""
    spirals(inner, outer; theta0=0., thickness=0.014, myelin=false, chi_I=-0.1, chi_A=-0.1, positions=[0, 0], repeats=[Inf, Inf], rotation=I(3)

Creates one or more [`Spiral`](@ref).
Spirals range from `inner` to `outer` radii starting at an angle of `theta0`.
Each wrap has a thickness of `thickness` micrometers.
Inner/outer cylinders can be added using respectively the `inner_cylinder` and `outer_cylinder` flags (default: only inner cylinder).

[Myelinated spirals](@ref Myelinated_annuli) can be created by setting the `myelin` to true.
All parameters can be either a single value or a vector of values.

The `positions`, `repeats`, and `rotation` control the annulus position and orientation and is explained in 
more detail in [Defining the geometry](@ref).
Additional keyword arguments are available to set generic obstruction settings as described in [`ObstructionProperties`](@ref).
"""
function spirals(args...; kwargs...)
    TransformObstruction(Spiral, args...; kwargs...)
end


function detect_collision(movement :: Movement{2}, spiral :: Spiral, previous :: Collision)
    # check collisions with cylinders first
    if spiral.inner_cylinder
        c_inner = detect_collision(movement, spiral.equivalent_annulus.inner, previous)
        if c_inner.index == 1
            # hit inside of inner cylinder
            return c_inner
        end
    else
        c_inner = empty_collision
    end
    if spiral.outer_cylinder
        c_outer = detect_collision(movement, spiral.equivalent_annulus.outer, previous)
        if c_outer.index == 0
            # hit inside of outer cylinder
            return c_outer
        end
    else
        c_outer = empty_collision
    end
    c_cylinders = c_inner.distance < c_outer.distance ? c_inner : c_outer


    toskip = previous === empty_collision ? -1 : previous.index
    for dim in 1:2
        if (
            (movement.origin[dim] > spiral.outer && movement.destination[dim] > spiral.outer) ||
            (movement.origin[dim] < -spiral.outer && movement.destination[dim] < -spiral.outer)
        )
            return empty_collision
        end
    end
    rsq(x) = x[1] * x[1] + x[2] * x[2]

    diff = movement.destination - movement.origin
    rsq_origin = rsq(movement.origin)
    rsq_destination = rsq(movement.destination)

    max_rsq = max(rsq_origin, rsq_destination)
    if max_rsq < (spiral.inner * spiral.inner)
        return empty_collision
    end

    a = rsq(diff)
    b = 2 * (diff ⋅ movement.origin)
    min_rsq_dist = - b / (2 * a)
    if min_rsq_dist < 0 || min_rsq_dist > 1
        if rsq_origin < rsq_destination
            min_rsq = rsq_origin
            min_rsq_pos = movement.origin
            min_rsq_dist = zero(Float)
        else
            min_rsq = rsq_destination
            min_rsq_pos = movement.destination
            min_rsq_dist = one(Float)
        end
    else
        min_rsq_pos = @. min_rsq_dist * movement.destination + (1 - min_rsq_dist) * movement.origin
        min_rsq = rsq(min_rsq_pos)
    end
    if min_rsq > (spiral.outer * spiral.outer)
        return empty_collision
    end

    slope = spiral.thickness / 2π
    theta_range = spiral.theta_end - spiral.theta0

    function get_theta(pos, rsq; round_down=true, ignore_toskip=true)
        theta = atan(pos[2], pos[1]) - spiral.theta0
        inner_radius = slope * theta + spiral.inner
        if ignore_toskip || (toskip == -1)
            nwrap = div(sqrt(rsq) - inner_radius, spiral.thickness, round_down ? RoundDown : RoundUp)
        else
            nwrap = div(sqrt(rsq) - inner_radius, spiral.thickness, RoundNearest)
            if !round_down && (toskip == 0)
                nwrap += 1
            elseif round_down && (toskip == 1)
                nwrap -= 1
            end
        end
        theta + nwrap * 2π
    end

    # y = line_slope * x + line_shift
    line_slope = diff[2] / diff[1]
    line_shift = -line_slope * movement.origin[1] + movement.origin[2]

    function froot(theta :: Float)
        (s, c) = sincos(theta + spiral.theta0)
        (s - line_slope * c) * (slope * theta + spiral.inner) - line_shift
    end
    theta_min_rsq = get_theta(min_rsq_pos, min_rsq; round_down=false, ignore_toskip=!iszero(min_rsq_dist))

    function get_solution(theta_sol, index)
        radius = slope * theta_sol + spiral.inner
        x = radius * cos(theta_sol + spiral.theta0)
        y = radius * sin(theta_sol + spiral.theta0)
        if iszero(diff[1])
            dist = (y - movement.origin[2]) / diff[2]
        else
            dist = (x - movement.origin[1]) / diff[1]
        end
        if c_cylinders.distance < dist
            return c_cylinders
        end
        return Collision(
            dist,
            iszero(index) ? SA[x, y, 0] : SA[-x, -y, 0],  # TODO: include spiralling in normal calculation
            ObstructionProperties(spiral),
            index=index
        )
    end

    function get_zero_solution()
        t1, t2 = get_theta(movement.origin, rsq_origin, ignore_toskip=false), get_theta(movement.destination, rsq_destination, ignore_toskip=true)
        if id(previous) == id(spiral)
            index = previous.index
            normal = previous.normal
        else
            index = t1 > t2 ? 1 : 0
            normal = get_solution(t1, index).normal
        end
        return Collision(zero(Float), normal, spiral.properties, previous.index)
    end

    if !iszero(min_rsq_dist) & (rsq_origin > spiral.inner * spiral.inner)
        # check downward trajectory
        theta_orig = get_theta(movement.origin, rsq_origin, ignore_toskip=false)
        if iszero(min_rsq)
            theta_min_rsq = mod(theta_orig, 2π) - 2π
        end
        theta_shift = mod((theta_min_rsq - theta_orig) + π, 2π) - π
        pos_shift = sign(theta_shift) > 0
        for theta1 in theta_orig:-2π:(theta_min_rsq - theta_shift - π)
            if abs(theta_shift) < 1e-8
                # particle heading straight for the centre
                if 0 < theta1 < theta_range
                    return get_solution(theta1, 0)
                end
                continue
            end
            (lower, upper) = pos_shift ? (theta1, theta1 + theta_shift) : (theta1 + theta_shift, theta1)
            if (lower > theta_range)
                continue
            end
            if (upper > theta_range)
                if sign(froot(lower)) == sign(froot(theta_range))
                    continue
                end
                upper = theta_range
            end
            if sign(froot(lower)) == sign(froot(upper))
                return get_zero_solution()
            end
            theta_sol = Roots.find_zero(froot, (lower, upper))
            if (theta_sol < 0)
                break
            end
            return get_solution(theta_sol, 0)
        end
    end

    if (rsq_destination > spiral.inner * spiral.inner)
        # check upward trajectory
        theta_dest = get_theta(movement.destination, rsq_destination, ignore_toskip=true)
        if iszero(min_rsq)
            theta_ref = mod(theta_dest, 2π) - 2π
        else
            theta_ref = get_theta(movement.origin, rsq_origin, ignore_toskip=false, round_down=false)
            if theta_ref > theta_range
                theta_ref -= div(theta_ref-theta_range, 2π, RoundDown) * 2π
            end
        end
        theta_shift = mod((theta_dest - theta_ref) + π, 2π) - π
        pos_shift = sign(theta_shift) > 0
        for theta1 in theta_ref:2π:(theta_dest - theta_shift + π)
            if abs(theta_shift) < 1e-8
                # particle heading straight out of centre
                if 0 < theta1 < theta_range
                    return get_solution(theta1, 1)
                end
                continue
            end
            (lower, upper) = pos_shift ? (theta1, theta1 + theta_shift) : (theta1 + theta_shift, theta1)
            if (upper < 0)
                continue
            end
            if (lower < 0)
                if sign(froot(upper)) == sign(froot(zero(Float)))
                    continue
                end
                lower = 0
            end
            if sign(froot(lower)) == sign(froot(upper))
                return get_zero_solution()
            end
            theta_sol = Roots.find_zero(froot, (lower, upper))
            if (theta_sol > theta_range)
                break
            end
            return get_solution(theta_sol, 1)
        end
    end
    t1, t2 = get_theta(movement.origin, rsq_origin, ignore_toskip=false), get_theta(movement.destination, rsq_destination, ignore_toskip=true)
    if abs(t1 - t2) > π
        return get_zero_solution()
    end
    return empty_collision
end

function lorentz_off_resonance(spiral::Spiral, position::SVector{2, Float}, b0_field::SVector{2, Float}, repeat_dist::SVector{2, Float}, radius::Float, nrepeats::SVector{2, Int})
    lorentz_off_resonance(spiral.equivalent_annulus, position, b0_field, repeat_dist, radius, nrepeats)
end

"""
    random_spirals(target_density; repeats, g_ratio=0.8, distribution=Distributions.Gamma, mean_radius=1., variance_radius=0.5, max_iter=1000, rotation=I(3))

Generate infinitely repeating box with non-overlapping spirals.

A rectangle with the size of `repeats` will be filled with spirals for a total surface density of `target_density`.
The spiral outer radii will be drawn from the selected `distribution` (if not set, a Gamma distribution is used with given `mean_radius` and `var_radius`).
An error is raised if no solution for non-overlapping annuli is found.
The inner radius with respect to the outer radius is set by the `g-ratio`.
Other spiral parameters (besides `inner`, `outer`, `positions`, and `repeats`) are identical as in `mr.spirals`.
"""
function random_spirals(target_density; repeats, g_ratio=0.8, distribution=nothing, mean_radius=1., variance_radius=0.5, max_iter=1000, kwargs...)
    (positions, outer) = random_positions_radii(repeats, target_density, 2; distribution=distribution, mean=mean_radius, variance=variance_radius, max_iter=max_iter)
    inner = g_ratio .* outer
    spirals(inner, outer; positions=positions, repeats=repeats, kwargs...)
end