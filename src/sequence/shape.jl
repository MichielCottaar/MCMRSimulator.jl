"""
    Shape(times, amplitudes; normalise=false)

Defines a 1D [`Shape`](@ref) for an RF pulse or gradient profile.
Times should be between 0 and 1 (inclusive).
Amplitudes should start and end at 0 and range between 0 and 1 (for RF pulses) and -1 and 1 (for gradients).
Errors are ranged if times/amplitudes are outside of these ranges (although normalise=true can be used to simply adjust the shape to match these constraints).
"""
struct Shape{N}
    times :: SVector{N, Float}
    amplitudes :: SVector{N, Float}
    function Shape(times, amplitudes; normalise=false)
        if !normalise
            @assert isapprox(minimum(times), 0., atol=1e-8)
            @assert isapprox(maximum(times), 1., atol=1e-8)
            @assert maximum(abs.(amplitudes)) == 1 || all(iszero.(amplitudes))
        end
        times = [times...]
        times .-= times[1]
        times ./= times[end]
        if !all(iszero.(amplitudes))
            amplitudes ./= maximum(abs.(amplitudes))
        end
        @assert length(times) == length(amplitudes)
        @assert all(times[2:end] .>= times[1:end-1])
        new{length(times)}(times, amplitudes)
    end
end

get_index(shape::Shape, time::Number) = time > 1 ? 0 : searchsortedlast(shape.times, time)
control_points(shape::Shape) = shape.times

Base.length(shape::Shape{N}) where {N} = N

"""
    amplitude(shape, time)

Returns the amplitude of the shape at the given time using linear interpolation.
"""
function amplitude(shape::Shape, time::Number)
    index = get_index(shape, time)
    if iszero(index)
        return zero(Float)
    end
    dt0 = time - shape.times[index]
    if iszero(dt0)
        return shape.amplitudes[index]
    end
    dt1 = shape.times[index + 1] - time
    return (shape.amplitudes[index] * dt1 + shape.amplitudes[index + 1] * dt0) / (dt1 + dt0)
end

"""
    amplitude_derivative(shape, time)

Returns the derivative of the shape at the specific time.
At control points the derivative is always assumed to be the one of the next block.
"""
function amplitude_derivative(shape::Shape{N}, time::Number) where {N}
    index = get_index(shape, time)
    if iszero(index) || (index == N)
        return zero(Float)
    end
    return (shape.amplitudes[index + 1] - shape.amplitudes[index]) / (shape.times[index + 1] - shape.times[index])
end

"""
    amplitude_integral_step(shape, index, t0, t1)

Integrate the amplitude from t0 to t1 assuming that both are within the single step (given by `index`).
This is a helper function for [`amplitude_integral`](@ref).
"""
function amplitude_integral_step(shape::Shape, index::Int, t0=nothing, t1=nothing)
    a0 = isnothing(t0) ? shape.amplitudes[index] : amplitude(shape, t0)
    a1 = isnothing(t1) ? shape.amplitudes[index] : amplitude(shape, t1)
    t0_real = isnothing(t0) ? shape.times[index] : t0
    t1_real = isnothing(t1) ? shape.times[index+1] : t1
    mean_amplitude = (a0 + a1) / 2
    return mean_amplitude * (t1_real - t0_real)
end

"""
    amplitude_integral(shape, t0, t1)

Integrate the amplitude from t0 to t1.
"""
function amplitude_integral(shape::Shape, t0::Number, t1::Number)
    @assert t1 >= t0
    if (t0 > 1) || (t1 < 0)
        return zero(Float)
    end
    if t0 < 0
        t0 = zero(Float)
    end
    if t1 > 1
        t1 = one(Float)
    end
    index0 = get_index(shape, t0)
    index1 = get_index(shape, t1)
    if index0 == index1
        return amplitude_integral_step(shape, index0, t0, t1)
    end
    total = amplitude_integral_step(shape, index0, t0, nothing)
    for index in index0+1:index1-1
        total += amplitude_integral_step(shape, index)
    end
    total += amplitude_integral_step(shape, index1, nothing, t1)
    return total
end

"""
    amplitude_integral(shape, t0, times)

Integrate the amplitude from t0 to all the times in `times` (all assumed to be larger than `t0` and strictly increasing).
"""
function amplitude_integral(shape::Shape, t0::T, times::AbstractVector{T}) where {T <: Number}
    totals = zeros(Float, length(times))
    current = zero(Float)
    prev_time = t0
    for (index, t) in enumerate(times)
        current += amplitude_integral(shape, prev_time, t)
        totals[index] = current
    end
    return totals
end


"""
    ConcreteShape(t0, t1, max_amplitude, shape)

A 1-dimensional amplitude profile extending from `t0` to `t1` with maximum amplitude of `max_amplitude`.
The amplitude profile is defined by a [`Shape`](@ref).
"""
struct ConcreteShape{N}
    t0::Float
    t1::Float
    max_amplitude::Float
    shift::Float
    shape::Shape{N}
    ConcreteShape(t0::Number, t1::Number, max_amplitude::Number, shift::Number, shape::Shape{N}) where {N} = new{N}(Float(t0), Float(t1), Float(max_amplitude), Float(shift), shape)
end

"""
    ConcreteShape(times, amplitudes)

Creates an amplitude profile by linearly interpolating the amplitudes between the given times.
Times and amplitudes should be vectors of the same length.
"""
function ConcreteShape(times::AbstractVector{<:Number}, amplitudes::AbstractVector{<:Number})
    if length(times) == 0
        return ConcreteShape([0, 1], [0, 0])
    end
    if length(times) <= 1
        error("ConcreteShape needs at least 2 control points")
    end
    t0 = times[1]
    t1 = times[end]
    max_amplitude = maximum(abs.(amplitudes))
    ConcreteShape(
        t0, t1, max_amplitude, 0,
        Shape(
            (times .- t0) ./ (t1 - t0),
            iszero(max_amplitude) ? amplitudes : (amplitudes ./ max_amplitude)
        )
    )
end

"""
    ConcreteShape()

Produces a ConcreteShape with only zero amplitude.
"""
ConcreteShape() = ConcreteShape(Float[], Float[])

convert_time(concrete::ConcreteShape, time::Number) = (time - concrete.t0) / (concrete.t1 - concrete.t0)
amplitude(concrete::ConcreteShape, time::Number) = amplitude(concrete.shape, convert_time(concrete, time)) * concrete.max_amplitude + concrete.shift
amplitude_derivative(concrete::ConcreteShape, time::Number) = amplitude_derivative(concrete.shape, convert_time(concrete, time)) * concrete.max_amplitude / (concrete.t1 - concrete.t0)
amplitude_integral(concrete::ConcreteShape, t0::Number, t1::Number) = (amplitude_integral(concrete.shape, convert_time(concrete, t0), convert_time(concrete, t1)) * concrete.max_amplitude + concrete.shift) * (concrete.t1 - concrete.t0)
amplitude_integral(concrete::ConcreteShape, t0::T, times::AbstractVector{T}) where {T<:Number} = amplitude_integral(concrete.shape, convert_time(concrete, t0), convert_time.(concrete, times)) .* (concrete.max_amplitude * (concrete.t1 - concrete.t0)) .+ (concrete.shift * (concrete.t1 - concrete.t0))
amplitude_integral(concrete::ConcreteShape, t0::T, times::AbstractRange{T}) where {T<:Number} = amplitude_integral(concrete.shape, convert_time(concrete, t0), convert_time(concrete, times)) .* (concrete.max_amplitude * (concrete.t1 - concrete.t0)) .+ (concrete.shift * (concrete.t1 - concrete.t0))
control_points(concrete::ConcreteShape) = control_points(concrete.shape) .* (concrete.t1 - concrete.t0) .+ concrete.t0

add_TR(concrete::ConcreteShape, TR::Number) = ConcreteShape(concrete.t0 + TR, concrete.t1 + TR, concrete.max_amplitude, concrete.shift, concrete.shape)

start_time(concrete::ConcreteShape) = concrete.t0
end_time(concrete::ConcreteShape) = concrete.t1


"""
    ShapePart(concrete_shape, t0, t1)

Represents a small part of a [`ConcreteShape`](@ref) between `t0` and `t1` during which the amplitude varies linearly.
This object is used to store the relevant part of the RF and gradient profiles within a single timestep (see [`SequencePart`](@ref)).
"""
struct ShapePart
    start :: Float
    final :: Float
    slope :: Float
end

function ShapePart(concrete_shape::ConcreteShape, t0::Number, t1::Number)
    tmean = (t0 + t1) / 2
    mean_value = amplitude(concrete_shape, tmean)
    slope = amplitude_derivative(concrete_shape, tmean)
    return ShapePart(
        (t0 - tmean) * slope + mean_value,
        (t1 - tmean) * slope + mean_value,
        slope * (t1 - t0),
    )
end

amplitude(part::ShapePart, time) = part.start + time * part.slope
amplitude(part::ShapePart, t1, t2) = amplitude(part, (t1 + t2) / 2)