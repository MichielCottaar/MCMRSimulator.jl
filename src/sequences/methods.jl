"""
Defines shared methods in the Sequence sub-module.
"""
module Methods

"""
    start_time(gradient/pulse)

Find the time that the MR gradient or pulse starts.
For instantaneous gradients/pulses this will also be [`end_time`](@ref).
"""
function start_time end

"""
    end_time(gradient/pulse)

Find the time that the MR gradient or pulse ends.
For instantaneous gradients/pulses this will also be [`start_time`](@ref).
"""
function end_time end

"""
    add_TR(shape/gradient/rf_pulse, TR)

Shifts the generic `Shape`, `MRGradients`, or `RFPulse` by a time `TR`.
"""
function add_TR end

end