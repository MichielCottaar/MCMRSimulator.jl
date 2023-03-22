"""
    MRGradients(times, amplitudes; origin=[0, 0, 0])
    MRGradients([(time0, amplitude0), (time1, amplitude1), ...]; origin=[0, 0, 0])

Defines a gradient profile with the gradients (unit: kHz/um) linearly interpolated between the given times (unit: ms).
Amplitudes can be a 3D vector or a single value. In the latter case the gradients are assumed to point in the z-direction (can be rotated using [`rotate_bvec`](@ref)).
The gradients are centered on given `origin` (unit: um).
They can be sampled using [`gradient`](@ref).

    MRGradients(Gx, Gy, Gz; origin=[0, 0, 0])

Builds the MR gradients out of 3 [`Shape`](@ref) objects each describing the gradient profile in 1 dimension.
"""
struct MRGradients
    shape::Shape{PosVector}
    origin::PosVector
end

function MRGradients()
    MRGradients(
        Shape{PosVector}(),
        zero(PosVector)
    )
end

function MRGradients(controls::AbstractVector{<:Tuple{<:Number, <:Any}}; origin=zero(PosVector))
    times = [c[1] for c in controls]
    amplitudes = [c[2] for c in controls]
    MRGradients(times, amplitudes; origin=origin)
end

function MRGradients(Gx::Shape, Gy::Shape, Gz::Shape; origin=zero(PosVector))
    gradients = [Gx, Gy, Gz]
    times = sort(unique(vcat(control_points.(gradients)...)))
    amplitudes = [PosVector([sample(G, t) for G in gradients]) for t in times]
    return MRGradients(times, amplitudes; origin=origin)
end

function MRGradients(times::AbstractVector{<:Number}, amplitudes::AbstractVector{<:Number}; origin=zero(PosVector))
    MRGradients(
        times,
        [mr.PosVector(0, 0, a) for a in amplitudes],
        origin=PosVector(origin)
    )
end

function MRGradients(times::AbstractVector{<:Number}, amplitudes::AbstractVector{<:AbstractVector{<:Number}}; origin=zero(PosVector))
    MRGradients(
        Shape(times, PosVector.(amplitudes)),
        PosVector(origin)
    )
end

add_TR(g::MRGradients, TR::Number) = MRGradients(add_TR(g.shape, TR), g.origin)

control_points(g::MRGradients) = control_points(g.shape)


"""
    gradient([position, ]sequence/grad, time)
    gradient([position, ]sequence/grad, t1, t2)

Get gradient strength generated by `grad` ([`MRGradients`](@ref) object) at a specific `time` or averaged over a period between `t1` and `t2`.
If no position is provided, the gradient is returned as a length-3 vector in units of mT/m.
If a position is provided, the gradient is returned at that position as a float quantifying the off-resonance field in units of microTesla.
"""
gradient(grad::MRGradients, time::Number) = sample(grad.shape, time)
gradient(grad::MRGradients, t1::Number, t2::Number) = sample(grad.shape, t1, t2)

function gradient(position::PosVector, grad::MRGradients, times...)
    g = gradient(grad, times...)
    return (g ⋅ (position - grad.origin))
end


"""
    rotate_bvec(gradients, bvec)

Rotates the gradients in `gradients` to align with `bvec`.
"""
function rotate_bvec(gradients::AbstractVector{<:Tuple{Real, Real}}, bvec)
    bvec = PosVector(bvec ./ norm(bvec))
    [(time, bvec * grad) for (time, grad) in gradients]
end

function rotate_bvec(gradients::MRGradients, bvec)
    rotation = get_rotation(bvec, 3)
    MRGradients(
        Shape(
            gradients.shape.times,
            [rotation * g for g in gradients.shape.amplitudes],
        ),
        gradients.origin
    )
end

