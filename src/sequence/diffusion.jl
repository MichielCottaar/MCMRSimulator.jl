const Optional{T} = Union{Nothing, T}

"""
    derive_qval_time(max_diffusion_time; bval=nothing, diffusion_time=nothing, qval=nothing)

Derives b-value, diffusion time, and q-value from each other assuming short pulses

Returns a tuple with the q-value and diffusion time.

The strength of the diffusion weighting can be defined in 3 ways:
- not setting anything or only setting the diffusion time: in this case no diffusion weighting is applied (b=0).
- setting only the b-value or q-value: in this case a `diffusion_time` of TE/2 is assumed.
- setting two out of three of the `b-value`, `diffusion_time`, and `qval`. The last one will be calculated.
If all three are defined an AssertionError is raised if they do not agree with each other.

Assumes equation of ``b = q^2 \\Delta``, where
- `bval` is the diffusion-weighted strength (b-value) in ``ms/um^2``.   
- `qval` is the gradient applied due to the diffusion-weighted gradients in the spin phase field (``rad/\\mu m``). 
    For a square pulse this is computed as ``\\gamma G \\delta``.
- `diffusion_time` (``\\Delta``) is the diffusion time defined for infinitely short pulses as the time between the diffusion-weighted gradients (in ``ms``).
"""
function derive_qval_time(
    bval :: Optional{Real},
    diffusion_time :: Optional{Real},
    qval :: Optional{Real},
    max_diffusion_time :: Real,
)
    if isnothing(bval)
        if isnothing(qval)
            qval = 0.
        end
        if isnothing(diffusion_time)
            diffusion_time = max_diffusion_time / 2.
        end
    elseif isnothing(qval)
        if isnothing(diffusion_time)
            diffusion_time = max_diffusion_time / 2.
        end
        qval = sqrt(bval / diffusion_time)
    elseif isnothing(diffusion_time)
        diffusion_time = bval / (qval * qval)
        if diffusion_time >= max_diffusion_time
            error("Implied diffusion time exceeds the maximal possible value.")
        end
    else
        @assert bval ≈ (qval * qval) * diffusion_time
    end
    (qval, diffusion_time)
end

derive_qval_time(
    max_diffusion_time :: Real;
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
) = derive_qval_time(bval, diffusion_time, qval, max_diffusion_time)


"""
    perfect_dwi(; TE=80., TR=2000., bval=nothing, diffusion_time=nothing, qval=nothing, orientation = [0, 0, 1.])

Creates a standard diffusion-weighted MRI sequence

This produces a pulsed gradient spin echo (PGSE) sequence with perfect RF pulses 
and infinitely short diffusion-weighted gradients.

The procedure for defining the diffusion-weighted gradient strength is defined by [`derive_qval_time`](@ref).
The `orientation` parameter only sets the gradient orientation, not the strength.
"""
function perfect_dwi(;
    TE=80.,
    TR=2000.,
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
    orientation=SVector{3, Float}([1., 0., 0.]),
)
    (qval, diffusion_time) = derive_qval_time(bval, diffusion_time, qval, TE)
    @assert diffusion_time < TE
    base_components = SequenceComponent[
        RFPulse(time=0., flip_angle=90., phase=-90.),
        RFPulse(time=TE/2., flip_angle=180.),
        Readout(time=TE),
    ]
    if !iszero(qval)
        qvec = PosVector(orientation .* (qval / norm(orientation)))
        append!(base_components, [
            InstantGradient(time=(TE - diffusion_time) / 2., qvec=qvec),
            InstantGradient(time=(TE + diffusion_time) / 2., qvec=qvec),
        ])
    end
    Sequence(pulses=base_components, TR=TR)
end

"""
    dwi(; TE=80., TR=2000., bval/qval/gradient_strength, diffusion_time=TE/2, gradient_duration, readout_time=0., scanner=Scanner(B0=3.), orientation = [1., 0, 0])

Creates a diffusion-weighted pulsed gradient spin echo MRI sequence

The timings of the RF pulses is set by `TE` and `TR`. 
The gradient timings will also be affected by `gradient_duration`, `diffusion_time`, `scanner` (which sets the ramp time) and `readout_time`:
- By default the gradient durations are set to the maximum value possible within the echo time (`TE`) keeping in mind the time needed for the MR readout (`readout_time`) and the time needed to ramp to the maximum gradient strength (set by the `scanner`).
- When `gradient_duration` is set to 0, the gradient pulses are assumed to be instanteneous (i.e., using [`InstantGradient`](@ref)). The time between these instant gradients can be set using `diffusion_time` (defaults to `TE`/2).

The strength of the diffusion gradients is set by one of `bval`, `qval`, or `gradient_strength`.
If this strength exceeds the maximum allowed for the `scanner` an AssertionError is raised.
The gradient orientation is set by `orientation`.
"""
function dwi(;
    TE=80.,
    TR=2000.,
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
    gradient_strength=nothing, # in mT/m
    gradient_duration=nothing,
    readout_time=0.,
    orientation=SVector{3, Float}([1., 0., 0.]),
    scanner=Scanner(B0=3.)
)
    if !isnothing(gradient_duration) && iszero(gradient_duration)
        return perfect_dwi(TE=TE, TR=TR, bval=bval, diffusion_time=diffusion_time, qval=qval, orientation=orientation)
    end

    grad_1D = dwi_gradients_1D(TE=TE, bval=bval, qval=qval, gradient_strength=gradient_strength, gradient_duration=gradient_duration, diffusion_time=diffusion_time, readout_time=readout_time, scanner=scanner)
    grad = rotate_bvec(grad_1D, orientation)
    pulses = SequenceComponent[
        RFPulse(time=0., flip_angle=90., phase=-90.),
        RFPulse(time=TE/2., flip_angle=180.),
        Readout(time=TE),
    ]
    Sequence(scanner=scanner, pulses=pulses, gradients=grad, TR=TR)
end


function dwi_gradients_1D(;
    TE,
    bval=nothing,
    qval=nothing,
    gradient_strength=nothing, # in mT/m
    gradient_duration=nothing, # gradient_duration defined as the duration where gradient strength equals the gradient_strength, the same as in the equation, ALL TIME ARGS ARE IN MILLISECOND.
    diffusion_time=nothing,
    readout_time=0.,
    scanner=Scanner(B0=3.)
)
    # determine ramp time, diffusion time and gradient duration and the timing of gradients
    t1 = 0.
    t2 = TE/2
    if isinf(scanner.gradient) || isinf(scanner.slew_rate) # Maybe implement a check for only inf gradient max in scanner()?
        ramp_time = 0.
    else
        ramp_time = scanner.gradient / scanner.slew_rate
    end

    if isnothing(diffusion_time)
        diffusion_time = TE/2
        if isnothing(gradient_duration)
            gradient_duration = (TE - readout_time)/2 - ramp_time
        else
            if gradient_duration == 0
                if isnothing(gradient_strength)
                    error("This case should have been dealt with in dwi, if you see this there's a problem")
                else
                    error("Can't have defined gradient strength when using instantaneous gradients")
                end
            elseif gradient_duration + diffusion_time + ramp_time + readout_time/2 <= TE
                t1, t2 = fit_time(gradient_duration=gradient_duration, diffusion_time=diffusion_time, ramp_time=ramp_time, readout_time=readout_time, TE=TE)
            else
                error("The timings specified cannot be fit within the given TE") # Since ramp time here is calculated from (max grad/ max slew rate), there may exist cases where the function will regard realistic sequences as unrealistic.
            end
        end
    else
        if isnothing(gradient_duration)
            if diffusion_time + ramp_time + readout_time/2 <= TE
                gradient_duration = min(TE - (diffusion_time + ramp_time + readout_time/2), diffusion_time, (TE - readout_time)/2 - ramp_time)
                t1, t2 = fit_time(gradient_duration=gradient_duration, diffusion_time=diffusion_time, ramp_time=ramp_time, readout_time=readout_time, TE=TE)
            else
                error("No time is left for applying the gradient")
            end
        elseif TE >= gradient_duration + diffusion_time + ramp_time + readout_time/2
            t1, t2 = fit_time(gradient_duration=gradient_duration, diffusion_time=diffusion_time, ramp_time=ramp_time, readout_time=readout_time, TE=TE)
        else
            error("TE is too short to fit everything in")
        end
    end
    
    if ramp_time > gradient_duration
        error("Gradient duration needs to be longer than or equal to the ramp time")
    end



    pulse_duration = gradient_duration
    if sum(isnothing.([bval, qval, gradient_strength])) != 2
        error("One and only one of the bval, qval, grdient_strength has to be defined")
    end
    
    if !isnothing(bval)
        gradient_strength = 1e3*sqrt(bval/(((diffusion_time - gradient_duration/3)*gradient_duration^2 + (ramp_time^3)/20 - (gradient_duration*ramp_time^2)/6)*(gyromagnetic_ratio^2)))
    elseif !isnothing(qval)
        gradient_strength = 1e3*qval/gradient_duration/gyromagnetic_ratio
    end
    @assert gradient_strength <= scanner.gradient "Requested gradient strength exceeds scanner limits"

    return [
        (t1, 0.),
        (t1 + ramp_time, gradient_strength),
        (t1 + pulse_duration, gradient_strength),
        (t1 + pulse_duration + ramp_time, 0.),
        (t2, 0.),
        (t2 + ramp_time, gradient_strength),
        (t2 + pulse_duration, gradient_strength),
        (t2 + pulse_duration + ramp_time, 0.),
    ]
end


""" 
    fit_time(gradient_duration, diffusion_time, ramp_time, readout_time, TE)

Arranges the timing of gradients given the duration of gradient, diffusion time, ramp time, readout time and TE. By default it will try to arrange to gradients symmetrically arround the 180 degree pulse.
If not feasible, it will make readout right after the rephasing gradient and calculate the dephasing gradient's timing based on the rephasing gradient.

This function is not supposed to be called individually! So no check for the arguments will happen, assuming all have been checked by dwi_gradients_1D and all will not be nothing.

"""
function fit_time(;
    gradient_duration, 
    diffusion_time, 
    ramp_time, 
    readout_time, 
    TE
) 
    if (diffusion_time + gradient_duration + ramp_time)/2 <= (TE - readout_time)/2 # Check the feasibility of symmetric arrangement
        t1 = TE/2 - ((diffusion_time + gradient_duration + ramp_time)/2)
        t2 = t1 + diffusion_time
    else
        t2 = TE - readout_time/2 - ramp_time - gradient_duration
        t1 = t2 - diffusion_time
    end
    return [t1, t2]
end
