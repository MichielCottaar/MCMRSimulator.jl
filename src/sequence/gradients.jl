"""
    create_gradients([(<start_time>, <gradient>), ...], interpolate=:step/:linear)

Models MRI gradients as a sequence of stepwise or linear functions.

Off-resonance fields are estimated by calling [`get_gradients`](@ref).
"""
abstract type MRGradients end

"""
    get_gradient([position, ]sequence/grad, time)
    get_gradient([position, ]sequence/grad, t1, t2)

Get gradient strength generated by `grad` ([`MRGradients`](@ref) object) at a specific `time` or averaged over a period between `t1` and `t2`.
If no position is provided, the gradient is returned as a length-3 vector.
If a position is provided, the gradient is returned at that position as a float quantifying the off-resonance field.
"""
function get_gradient(position::AbstractVector{<:Real}, grad::MRGradients, time::Real)
    get_gradient(PosVector(position), grad, Float(time))
end

function get_gradient(position::PosVector, grad::MRGradients, time::Float)
    g = get_gradient(grad, time)
    return g ⋅ (position - grad.origin)
end

function get_gradient(position::AbstractVector{<:Real}, grad::MRGradients, t1::Real, t2::Real)
    get_gradient(PosVector(position), grad, Float(t1), Float(t2))
end

function get_gradient(position::PosVector, grad::MRGradients, t1::Float, t2::Float)
    g = get_gradient(grad, t1, t2)
    return g ⋅ (position - grad.origin)
end


"""
    create_gradients([(<start_time>, <gradient>), ...], interpolate=:step)

Models MRI gradients as a sequence of stepwise functions.
"""
struct StepWiseGradients <: MRGradients
    times::AbstractVector{Float}
    gradients::AbstractVector{PosVector}
    origin::PosVector
    function StepWiseGradients(steps, origin=zero(PosVector); TR)
        times = Float[Float(s[1]) for s in steps]
        pushfirst!(times, zero(Float))
        push!(times, Float(TR))
        @assert issorted(times)

        gradients = PosVector[PosVector(s[2]) for s in steps]
        pushfirst!(gradients, zero(PosVector))
        push!(gradients, zero(PosVector))

        new(times, gradients, PosVector(origin))
    end
end

function get_gradient(grad::StepWiseGradients, time::Real)
    if time == grad.times[end]
        return grad.gradients[end]
    end
    if time == grad.times[1]
        return grad.gradients[1]
    end

    index = searchsortedfirst(grad.times, time) - 1
    return grad.gradients[index]
end

function get_gradient(grad::StepWiseGradients, t1::Real, t2::Real)
    if t2 < t1
        # we have to simulate over a complete TR
        g1 = get_gradient(grad, t1, grad.times[end])
        g2 = get_gradient(grad, zero(Float), t2)
        dt1 = grad.times[end] - t1
        dt_total = dt1 + t2
        return @. (g1 * dt1 + g2 * dt2) / dt_total
    end
    i1 = searchsortedfirst(grad.times, t1) - 1
    i2 = searchsortedfirst(grad.times, t2) - 1
    if i1 == i2
        return grad.gradients[i1]
    end
    tnext = t1
    total = zeros(Float, 3)
    for index in i1:i2
        tprev = tnext
        tnext = index == i2 ? t2 : grad.times[index + 1]
        if tnext != tprev
            total = total .+ grad.gradients[index] .* (tnext - tprev)
        end
    end
    return PosVector(total) ./ (t2 - t1)
end


"""
    create_gradients([(<start_time>, <gradient>), ...], interpolate=:linear)

Models MRI gradients as a sequence of linear functions.
"""
struct LinearGradients <: MRGradients
    times::AbstractVector{Float}
    gradients::AbstractVector{PosVector}
    origin::PosVector
    function LinearGradients(steps, origin=zero(PosVector); TR)
        @assert issorted(steps, by=s->s[1])
        times = [Float(s[1]) for s in steps]
        pushfirst!(times, zero(Float))
        push!(times, Float(TR))
        gradients = [PosVector(s[2]) for s in steps]
        pushfirst!(gradients, zero(PosVector))
        push!(gradients, zero(PosVector))
        new(times, gradients, PosVector(origin))
    end
end

function get_gradient(grad::LinearGradients, time::Real)
    if time == grad.times[end]
        return grad.gradients[end]
    end
    if time == grad.times[1]
        return grad.gradients[1]
    end

    index = searchsortedfirst(grad.times, time) - 1
    t1 = grad.times[index]
    t2 = grad.times[index + 1]
    dt1 = (time - t1) / (t2 - t1)
    dt2 = 1 - dt1
    g1 = grad.gradients[index]
    g2 = grad.gradients[index + 1]
    return @. g2 * dt1 + g1 * dt2
end

function get_gradient(grad::LinearGradients, t1::Real, t2::Real)
    if t2 < t1
        # we have to simulate over a complete TR
        g1 = get_gradient(grad, t1, grad.times[end])
        g2 = get_gradient(grad, zero(Float), t2)
        dt1 = grad.times[end] - t1
        dt_total = dt1 + t2
        return @. (g1 * dt1 + g2 * t2) / dt_total
    end
    i1 = searchsortedfirst(grad.times, t1) - 1
    i2 = searchsortedfirst(grad.times, t2) - 1
    if i1 == i2
        return get_gradient(grad, (t1 + t2)/2)
    end
    tnext = t1
    total = zeros(Float, 3)
    for index in i1:i2
        tprev = tnext
        tnext = index == i2 ? t2 : grad.times[index + 1]
        if tnext != tprev
            total = total .+ get_gradient(grad, (tprev + tnext)/2) .* (tnext - tprev)
        end
    end
    return PosVector(total) ./ (t2 - t1)
end

"""
    create_gradients([(<start_time>, <gradient>), ...], TR; interpolate=:linear)

Models MRI gradients as a sequence of linear functions.
`interpolate` can be set to :step to model the gradients as stepwise functions.
"""
function create_gradients(steps::AbstractVector, TR::Real; origin=zero(PosVector), interpolate=:linear)
    if interpolate == :linear
        return LinearGradients(steps, origin; TR=TR)
    elseif interpolate == :step
        return StepWiseGradients(steps, origin; TR=TR)
    else
        error("interpolate should be :step or :linear, not $interpolate")
    end
end


"""
    rotate_bvec(gradients, bvec; orig_bvec=nothing)

Rotates the gradients in `gradients` to align with `bvec`.
If not provided, the original gradient orientation is estimated as the first eigenvector of the b-tensor.
"""
function rotate_bvec(gradients::AbstractVector{<:Tuple{Real, Real}}, bvec; orig_bvec=nothing)
    bvec = PosVector(bvec / norm(bvec))
    [(time, bvec * grad) for (time, grad) in gradients]
end

