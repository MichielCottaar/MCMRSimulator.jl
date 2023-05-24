"""
    dwssfp(TR; TE=<TR>, scanner=<3T scanner>, excitation_pulse=Instant, excitation_time=<half pulse duration>, flip_angle=90, gradient_strength=0, gradient_duration=0, gradient_time=nothing)

Creates a unit of a diffusion weight steady state free precession sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- a diffusion weighting gradient
- a readout `TE` ms after the excitation.
"""

function dwssfp(TR; 
                TE=nothing, 
                scanner=Scanner(B0=3.), 
                excitation_pulse=InstantRFPulse(flip_angle=90), 
                excitation_time=nothing, 
                flip_angle=90, 
                gradient_strength=0, 
                gradient_duration=0, 
                gradient_delay=nothing)



end