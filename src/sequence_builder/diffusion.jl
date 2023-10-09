module Diffusion
import ...Scanners: max_gradient, max_slew_rate, Scanner
import ...Sequences: MRGradients, rotate_bvec, InstantGradient
import ...Methods: get_rotation
import ..BuildingBlocks: isempty_block, duration, BuildingBlock

"""
    add_linear_diffusion_weighting(blocks, replace1, replace2, refocus=true, bval/qval/gradient_strength, diffusion_time=max, gradient_duration=max, scanner=3T, orientation=:x)

Replaces two empty [`BuildingBlock`](@ref)-like objects with diffusion-weighting gradients.
Both diffusion-weighted gradients are identical trapezia. The second one will be flipped if `refocus` is set to false.

The other flags control the timings and size of the diffusion weighting.
The timings are determined by:
- `diffusion_time`: defaults to time between the start of the two blocks
- `gradient_duration`: defaults to the largest value possible. Setting the `gradient_duration` to zero will cause [`InstantGradient`] objects to be used.
- `ramp_time`: defaults to the maximum gradient / maximum slew rate of the scanner.

The gradient strength is set by `bval` (ms/um^2), `qval` (rad/um), or `gradient_strength` (kHz/um). At least one of these should be provided.

By default, the timings of the existing [`BuildingBlock`](@ref) objects are respected.
If `max_bval` or `max_qval` are set these timings are adjusted to the minimum duration required to reach that b-value or q-value.
"""
function add_linear_diffusion_weighting(
    blocks, replace1, replace2;
    refocus=true, kwargs...
)
    @assert isempty_block(blocks[replace1])
    @assert isempty_block(blocks[replace2])
    @assert replace2 > replace1
    @assert (replace2 - replace1) > 1 || !refocus

    time_in_between = sum(duration.(blocks[replace1+1:replace2-1]))
    t1 = duration(blocks[replace1])
    t2 = duration(blocks[replace2])
    (trap1, trap2) = get_diffusion_trapeziums(t1, time_in_between, t2; kwargs...)
    if !refocus
        trap2[2] = MRGradients(
            Shape(
                trap2[2].shape.times,
                -trap2[2].shape.amplitudes,
            ),
            trap2[2].shape.PosVector
        )
    end

    new_blocks = copy(blocks)
    new_blocks[replace1] = trap1
    new_blocks[replace2] = trap2
    return new_blocks
end


function get_diffusion_trapeziums(
    duration_trap1, duration_wait, duration_trap2;
    bval=nothing, qval=nothing, gradient_strength=nothing, scanner=Scanner(B0=3.),
    diffusion_time=nothing, gradient_duration=nothing, ramp_time=nothing,
    orientation=:x,
)
    if isnothing(ramp_time)
        if isinf(max_gradient(scanner)) || isinf(max_slew_rate(scanner)) # Maybe implement a check for only inf gradient max in scanner()?
            ramp_time = 0.
        else
            ramp_time = max_gradient(scanner) / max_slew_rate(scanner)
        end
    end

    if isnothing(diffusion_time)
        diffusion_time = max(duration_trap1, duration_trap2) + duration_wait
    end

    if isnothing(gradient_duration)
        gradient_duration = min(
            duration_trap1,
            diffusion_time - duration_wait,
            duration_trap2,
            duration_trap2 - (diffusion_time - duration_trap1 - duration_wait)
        ) - ramp_time
        if gradient_duration < 0
            error("No time is left for applying the gradient")
        end
    elseif iszero(gradient_duration)
        ramp_time = 0.
    end

    t1, t2 = fit_time(
        duration_trap1, duration_wait, duration_trap2,
        gradient_duration, diffusion_time, ramp_time,
    )
    t1_after = duration_trap1 - gradient_duration - ramp_time - t1
    t2_after = duration_trap2 - gradient_duration - ramp_time - t2
    
    if ramp_time > gradient_duration
        error("Gradient duration needs to be longer than or equal to the ramp time")
    end

    if sum(isnothing.([bval, qval, gradient_strength])) != 2
        error("One and only one of the bval, qval, gradient_strength should be defined")
    end

    if iszero(gradient_duration)
        if !isnothing(gradient_strength)
            error("Cannot directly set `gradient_strength` if `gradient_duration` has been set to zero")
        elseif !isnothing(bval)
            qval = sqrt(bval / diffusion_time)
        end
        grad = InstantGradient(qvec=qval .* get_rotation(orientation, 1)[:, 1])
        return (
            [
                t1,
                grad,
                duration_trap1 - t1,
            ],
            [
                t2,
                grad,
                duration_trap2 - t2,
            ]
        )
    end
    
    if !isnothing(bval)
        gradient_strength = sqrt(bval/((diffusion_time - gradient_duration/3)*gradient_duration^2 + (ramp_time^3)/20 - (gradient_duration*ramp_time^2)/6)) / 2π
    elseif !isnothing(qval)
        gradient_strength = qval/gradient_duration / 2π
    end
    @assert gradient_strength <= max_gradient(scanner) "Requested gradient strength exceeds scanner limits"

    trapezium = rotate_bvec([
        (0, 0.), 
        (ramp_time, gradient_strength),
        (gradient_duration, gradient_strength),
        (gradient_duration + ramp_time, 0.), 
    ], orientation)
    return (
        [
            t1,
            trapezium,
            t1_after
        ],
        [
            t2,
            trapezium,
            t2_after
        ],
    )
end


""" 
    fit_time(duration1, duration_wait, duration2, gradient_duration, diffusion_time, ramp_time)

Arranges the timing of gradients given the duration of gradient, diffusion time, ramp time, readout time and TE. By default it will try to arrange to gradients symmetrically arround the 180 degree pulse.
If not feasible, it will make readout right after the rephasing gradient and calculate the dephasing gradient's timing based on the rephasing gradient.

This function is not supposed to be called individually! So no check for the arguments will happen, assuming all have been checked by `get_diffusion_trapeziums` and all will not be nothing.
"""
function fit_time(
    duration1, duration_wait, duration2,
    gradient_duration, diffusion_time, ramp_time, 
) 
    total_duration = diffusion_time + gradient_duration + ramp_time
    if total_duration/2 <= min(duration1, duration2) + duration_wait/2 # Check the feasibility of symmetric arrangement
        t1 = duration1 - (total_duration - duration_wait)/2
        t2 = (total_duration - duration_wait)/2 - gradient_duration - ramp_time
    elseif duration1 < duration2
        t1 = 0.
        t2 = diffusion_time - (duration1 + duration_wait)
    else
        t2 = 0.
        t1 = duration1 - (diffusion_time - duration_wait)
    end
    return (t1, t2)
end
end