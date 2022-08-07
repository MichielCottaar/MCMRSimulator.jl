const Optional{T} = Union{Nothing, T}

"""
Derives b-value, diffusion time, and q-value from each other assuming short pulses

Returns a tuple with the q-value and diffusion time.

The strength of the diffusion weighting can be defined in 3 ways:
- not setting anything or only setting the diffusion time: in this case no diffusion weighting is applied (b=0).
- setting only the b-value or q-value: in this case a `diffusion_time` of TE/2 is assumed.
- setting two out of three of the `b-value`, `diffusion_time`, and `qval`. The last one will be calculated.
If all three are defined an AssertionError is raised if they do not agree with each other.

Assumes equation of ``b = \\q^2 \\Delta``, where
- `bval` is the diffusion-weighted strength (b-value) in ``s/mm^2``.
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
Creates a standard diffusion-weighted MRI sequence

This produces a pulsed gradient spin echo (PGSE) sequence with perfect RF pulses 
and infinitely short diffusion-weighted gradients.

The procedure for defining the diffusion-weighted gradient strength is defined by `derive_qval_time`[@ref].
The `orientation` parameter only sets the gradient orientation, not the strength.
"""
function perfect_dwi(;
    TE=80.,
    TR=2000.,
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
    orientation=SVector{3, Float64}([0., 0., 1.]),
)
    (qval, diffusion_time) = derive_qval_time(bval, diffusion_time, qval, TE)
    @assert diffusion_time < TE
    base_components = SequenceComponent[
        RFPulse(time=0., flip_angle=90., phase=-90.),
        RFPulse(time=TE/2., flip_angle=180.),
        Readout(time=TE),
    ]
    if !iszero(qval)
        qvec = orientation .* (qval / norm(orientation))
        append!(base_components, [
            InstantGradient(time=(TE - diffusion_time) / 2., qvec=qvec),
            InstantGradient(time=(TE + diffusion_time) / 2., qvec=qvec),
        ])
    end
    Sequence(base_components, TR)
end