"""
    dwssfp(TR; TE=<TR>, scanner=<3T scanner>, excitation_pulse=Instant, excitation_time=<half pulse duration>, flip_angle=90, gradient_strength=0, gradient_duration=0, gradient_time=nothing)

Creates a unit of a diffusion weight steady state free precession sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- a diffusion weighting gradient (by default, the gradient will be centred at (TE-excitation_time)/2)
- a readout `TE` ms after the excitation.
"""

function dwssfp(TR; 
                TE=nothing, 
                scanner=Scanner(B0=3.), 
                excitation_pulse=InstantRFPulse(flip_angle=90), 
                excitation_time=0, 
                flip_angle=90, 
                gradient_strength=nothing, # Only used when gradient_duration != 0 and qval == nothing
                gradient_duration=0, 
                qval = nothing,
                gradient_delay=0.1, # have a small delay to separate RF pulse and gradient by default
                gradient_orientation=SVector{3, Float}([1., 0., 0.]))
    
    if isnothing(TE)
        TE = TR
    end

    if excitation_time == 0
        excitation_pulse = InstantRFPulse(flip_angle = flip_angle)
    else
        excitation_pulse = constant_pulse(0, excitation_time, flip_angle)
    end

    if gradient_duration == 0
        if qval == 0 
            gradient = []
        else
            gradient = InstantGradient(qvec=qval .* get_rotation(gradient_orientation, 1)[:, 1] / 2π) # ? ask Michiel
        end

        @assert gradient_delay < TE "gradient delay too long"
        ramp_time = 0
    else
        if isinf(max_gradient(scanner)) || isinf(max_slew_rate(scanner)) 
            ramp_time = 0.
        else
            ramp_time = max_gradient(scanner) / max_slew_rate(scanner)
        end
        
        if isnothing(qval) && isnothing(gradient_strength) 
            error("Either q value (preferred) or gradient strength needs to be specified for a finite gradient duration")
        elseif !isnothing(qval) && !isnothing(gradient_strength) 
            @assert gradient_strength == qval/gradient_duration / 2π   "gradient_strength and qval mismatch!"
        elseif !isnothing(qval) && isnothing(gradient_strength)
            gradient_strength = qval/gradient_duration / 2π 
        end
        @assert gradient_strength <= max_gradient(scanner) "Requested gradient strength exceeds scanner limits"
        gradient = rotate_bvec([
                (0, 0.), 
                (ramp_time, gradient_strength),
                (gradient_duration, gradient_strength),
                (gradient_duration + ramp_time, 0.), 
            ], gradient_orientation)

#        @assert gradient_delay > excitation_time
        if gradient_delay + gradient_duration + ramp_time + excitation_time > TE
            error("Can't fit gradient with the specified delay inside TE")
        end
    end
    # display(gradient)
    define_sequence(scanner, TR) do 
        [
            excitation_pulse,
            gradient_delay,
            gradient,
            TE - (gradient_delay + excitation_time + gradient_duration + ramp_time),
            Readout()
        ]
    end

end