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
    equivalent_annulus :: Annulus
    properties :: ObstructionProperties
end

function Spiral(inner::Float, outer::Float, theta0::Float, thickness::Float; kwargs...)
    @assert outer > inner
    closed = (outer - inner) > thickness
    theta_end = theta0 + Float(2π) * (outer - inner) / thickness
    Spiral(inner, outer, theta0, theta_end, thickness, closed, Annulus(inner, outer), ObstructionProperties(kwargs...))
end

Spiral(inner, outer; theta0=0., thickness=0.014, kwargs...) = Spiral(Float(inner), Float(outer), Float(theta0), Float(thickness); kwargs...)
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


function spirals(args...; kwargs...)
    TransformObstruction(Spiral, args...; kwargs...)
end


function detect_collision(movement :: Movement{2}, spiral :: Spiral, previous :: Collision)
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

    function get_theta(pos, rsq; round_down=true, ignore_toskip=true)
        theta = atan(pos[2], pos[1]) - spiral.theta0
        inner_radius = slope * theta + spiral.inner
        if ignore_toskip || (toskip == -1)
            nwrap = div(sqrt(rsq) - inner_radius, spiral.thickness, round_down ? RoundDown : RoundUp)
        else
            nwrap = div(sqrt(rsq) - inner_radius, spiral.thickness, RoundNearest)
            if !round_down && (toskip == 0)
                nwrap += 1
            elseif !round_down && (toskip == 1)
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
        dist = (x - movement.origin[1]) / diff[1]
        return Collision(
            dist,
            SA[-y, x, 0],  # TODO: include spiralling in normal calculation
            ObstructionProperties(spiral),
            index=index
        )
    end

    theta_range = spiral.theta_end - spiral.theta0

    if !iszero(min_rsq_dist) & (rsq_origin > spiral.inner * spiral.inner)
        # check downward trajectory
        theta_orig = get_theta(movement.origin, rsq_origin, ignore_toskip=false)
        theta_shift = mod((theta_min_rsq - theta_orig) + π, 2π) - π
        pos_shift = sign(theta_shift) > 0
        for theta1 in theta_orig:-2π:theta_min_rsq
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
            theta_sol = Roots.find_zero(froot, (lower, upper))
            @assert (theta_sol > theta_min_rsq) 
            if (theta_sol < 0)
                break
            end
            return get_solution(theta_sol, 0)
        end
    end

    if !isone(min_rsq_dist) & (rsq_destination > spiral.inner * spiral.inner)
        # check upward trajectory
        theta_dest = get_theta(movement.destination, rsq_destination, ignore_toskip=true)
        theta_shift = mod((theta_dest - theta_min_rsq) + π, 2π) - π
        pos_shift = sign(theta_shift) > 0
        for theta1 in theta_min_rsq:2π:theta_dest
            (lower, upper) = pos_shift ? (theta1, theta1 + theta_shift) : (theta1 + theta_shift, theta1)
            if (upper < 0)
                continue
            end
            if (lower < 0)
                if sign(froot(upper)) == sign(froot(0))
                    continue
                end
                lower = 0
            end
            theta_sol = Roots.find_zero(froot, (lower, upper))
            @assert (theta_sol > theta_min_rsq) 
            if (theta_sol > theta_range)
                break
            end
            return get_solution(theta_sol, 1)
        end
    end
    return empty_collision
end