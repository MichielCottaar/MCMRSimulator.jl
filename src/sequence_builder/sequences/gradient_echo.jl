"""
    gradient_echo(TE; TR=<TE>, scanner=<3T scanner>, excitation_pulse=Instant, excitation_time=<half pulse duration>)

Creates a gradient echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- a readout `TE` ms after the excitation.
"""
function gradient_echo(
    TE;
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90),
    excitation_time=nothing
)
    excitation_time = isnothing(excitation_time) ? duration(excitation_pulse) / 2 : excitation_time
    use_scanner(scanner) do 
        [
            excitation_pulse,
            TE - duration(excitation_pulse) + excitation_time,
            Readout(),
            TR - TE - excitation_time
        ]
    end
end