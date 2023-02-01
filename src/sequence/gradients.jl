"""
    MRGradients(Gx, Gy, Gz, origin)
    MRGradients(times, amplitudes; origin=[0, 0, 0])
    MRGradients([(time0, amplitude0), (time1, amplitude1), ...]; origin=[0, 0, 0])
    MRGradients()

Defines a gradient profile with the gradients (unit: kHz/um) linearly interpolated between the given times (unit: ms).
The gradients are centered on given `origin` (unit: um).
They can be sampled using `gradient`.
"""
struct MRGradients{Nx, Ny, Nz}
    Gx::ConcreteShape{Nx}
    Gy::ConcreteShape{Ny}
    Gz::ConcreteShape{Nz}
    origin::PosVector
end

function MRGradients()
    MRGradients(
        ConcreteShape(),
        ConcreteShape(),
        ConcreteShape(),
        zero(PosVector)
    )
end

function MRGradients(controls::AbstractVector{<:Tuple{<:Number, <:Any}}; origin=zero(PosVector))
    times = [c[1] for c in controls]
    amplitudes = [c[2] for c in controls]
    MRGradients(times, amplitudes; origin=origin)
end

function MRGradients(times::AbstractVector{<:Number}, amplitudes::AbstractVector{<:Number}; origin=zero(PosVector))
    MRGradients(
        ConcreteShape(times, amplitudes),
        ConcreteShape(),
        ConcreteShape(),
        PosVector(origin)
    )
end

function MRGradients(times::AbstractVector{<:Number}, amplitudes::AbstractVector{<:AbstractVector{<:Number}}; origin=zero(PosVector))
    @assert all([length(a) == 3 for a in amplitudes])
    shapes = [ConcreteShape(times, [a[index] for a in amplitudes]) for index in 1:3]
    MRGradients(
        shapes...,
        PosVector(origin)
    )
end

add_TR(g::MRGradients, TR::Number) = MRGradients(add_TR(g.Gx, Tr),add_TR(g.Gy, TR), add_TR(g.Gz, TR), g.origin)
control_points(g::MRGradients) = [
    control_points(g.Gx)...,
    control_points(g.Gy)...,
    control_points(g.Gz)...,
]


"""
    gradient([position, ]sequence/grad, time)
    gradient([position, ]sequence/grad, t1, t2)

Get gradient strength generated by `grad` ([`MRGradients`](@ref) object) at a specific `time` or averaged over a period between `t1` and `t2`.
If no position is provided, the gradient is returned as a length-3 vector in units of mT/m.
If a position is provided, the gradient is returned at that position as a float quantifying the off-resonance field in units of microTesla.
"""
function gradient(grad::MRGradients, time::Number)
    PosVector(
        amplitude(grad.Gx, time),
        amplitude(grad.Gy, time),
        amplitude(grad.Gz, time),
    )
end

function gradient(grad::MRGradients, t1::Number, t2::Number)
    PosVector(
        amplitude_integral(grad.Gx, t1, t2),
        amplitude_integral(grad.Gy, t1, t2),
        amplitude_integral(grad.Gz, t1, t2),
    ) / (t2 - t1)
end

function gradient(position::AbstractVector{<:Number}, grad::MRGradients, time::Number)
    gradient(PosVector(position), grad, Float(time))
end

function gradient(position::PosVector, grad::MRGradients, time::Float)
    g = gradient(grad, time)
    return (g ⋅ (position - grad.origin))
end

function gradient(position::AbstractVector{<:Number}, grad::MRGradients, t1::Number, t2::Number)
    gradient(PosVector(position), grad, Float(t1), Float(t2))
end

function gradient(position::PosVector, grad::MRGradients, t1::Float, t2::Float)
    g = gradient(grad, t1, t2)
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

