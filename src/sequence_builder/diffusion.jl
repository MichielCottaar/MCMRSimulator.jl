module Diffusion
import StaticArrays: SVector
import LinearAlgebra: norm
import ...Scanners: max_gradient, max_slew_rate, Scanner
import ...Sequences: MRGradients, rotate_bvec, InstantGradient
import ...Methods: get_rotation
import ..BuildingBlocks: isempty_block, duration, BuildingBlock

"""
    trapezium_gradient(; qval, total_duration, δ, gradient_strength, ramp_time, scanner, origin, apply_bvec, orientation)

Defines a MRI gradient shaped like a trapezium.

The resulting trapezium will obey the following properties:
- total_duration = δ + ramp_time
- qval = δ * gradient_strength
- gradient_strength < max_gradient(scanner)
- gradient_strength / ramp_time < max_slew_rate(scanner)

If the timings (`total_duration` or `δ`) are set, but the `qval` is not, then the `qval` will be the maximum possible within the gradient duration.
If the `qval` is set, but the timings are not, then the gradient duration is minimised.
In both cases the `gradient_strength` will be the maximum allowed by the scanner (if not explicitly set).

If the `ramp_time` is not explicitly set, it will be set to the slew rate over the gradient strength.

An [`InstantGradient`](@ref) is returned if `δ` or `total_duration` is set to zero (rather than an [`MRGradients`](@ref))
"""
function trapezium_gradient(; qval=nothing, total_duration=nothing, δ=nothing, gradient_strength=nothing, ramp_time=nothing, scanner=Scanner(B0=3.), origin=zero(SVector{3, Float64}), apply_bvec=false, orientation=[1., 0., 0.])
    orientation = orientation ./ norm(orientation)
    if apply_bvec
        maxima = (max_gradient(scanner), max_slew_rate(scanner))
    else
        maxima = (max_gradient(scanner) / maximum(abs.(orientation)), max_slew_rate(scanner) / maximum(abs.(orientation)))
    end

    if total_duration == 0 || δ == 0.
        if isnothing(qval)
            error("Need to set qval for trapezium gradient when having a gradient duration of 0")
        end
        if isnothing(total_duration) || iszero(total_duration)
            return InstantGradient(qvec=qval .* orientation, origin=origin, apply_bvec=apply_bvec)
        else
            return [total_duration / 2, InstantGradient(qvec=qval .* orientation, origin=origin, apply_bvec=apply_bvec), total_duration/2]
        end
    end
    if isinf(maxima[2])
        min_ramp_time = 0.
    else
        min_ramp_time = maxima[1] / maxima[2]
    end

    # determine timings
    if isnothing(total_duration)
        if isnothing(δ)
            if isnothing(qval)
                error("Need to set at least one of qval, total_duration, or δ when calling `get_trapezium`")
            end
            if isnothing(gradient_strength)
                gradient_strength = maxima[1]
                if isinf(gradient_strength)
                    return InstantGradient(qval .* orientation, origin, apply_bvec=apply_bvec)
                end
            end
            δ = qval / gradient_strength
        end
        if isnothing(ramp_time)
            ramp_time = min_ramp_time
        end
        total_duration = δ + ramp_time
    else
        if isnothing(δ)
            if isnothing(ramp_time)
                ramp_time = min_ramp_time
            end
            δ = total_duration - ramp_time
        else
            if isnothing(ramp_time)
                ramp_time = total_duration - δ
            end
        end
    end
    if δ < ramp_time
        ramp_time = total_duration / 2
        δ = total_duration / 2
    end


    # determine gradient strength
    if isnothing(qval)
        if isnothing(gradient_strength)
            if isinf(maxima[2])
                gradient_strength = maxima[1]
            else
                gradient_strength = min(maxima[1], ramp_time * maxima[2])
            end
            if isinf(gradient_strength)
                error("Trying to define a trapezium purely on its timings, however the q-value can be infintely largely because scanner does not have a max_gradient.")
            end
        end
        qval = gradient_strength * δ
    elseif isnothing(gradient_strength)
        gradient_strength = qval / δ
    end
    @assert ramp_time + δ ≈ total_duration
    @assert qval ≈ gradient_strength * δ
    @assert gradient_strength <= maxima[1] * 1.00001
    @assert iszero(ramp_time) || (gradient_strength / ramp_time <= maxima[2] * 1.00001)
    return MRGradients([0, ramp_time, δ, δ + ramp_time], [a * orientation for a in [0., gradient_strength, gradient_strength, 0.]], origin=origin, apply_bvec=apply_bvec)
end

"""
    gen_crusher(qval=<maximum>, duration=<minimum>, scanner=<scanner with infinite gradients>)

Generate a crusher gradient with the user-defined q-value (rad/um) and duration (ms).
If not provided, the q-value will be as large as possible, while the duration will be as small as possible given the constraints of the scanner (might be 0 for infinitely strong scanners).
For scanners with infinitely strong gradients (default) or durations of 0 ms, the q-value is set to 1 rad/um.
"""
function gen_crusher(; qval=nothing, duration=nothing, scanner=Scanner())
    if isnothing(qval) && (isinf(max_gradient(scanner)) || duration == 0)
        qval = 1.
    end
    return trapezium_gradient(total_duration=duration, scanner=scanner, qval=qval, orientation=[1, 1, 1], apply_bvec=false)
end

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

If the gradient orientation is not set during construction, it can be later applied using [`rotate_bvec`](@ref).
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
            trap2[2].shape.PosVector,
            true
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
        grad = InstantGradient(qvec=qval .* get_rotation(orientation, 1)[:, 1], apply_bvec=true)
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

    trapezium = rotate_bvec(MRGradients([
        (0, 0.), 
        (ramp_time, gradient_strength),
        (gradient_duration, gradient_strength),
        (gradient_duration + ramp_time, 0.), 
    ], apply_bvec=true), orientation)
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