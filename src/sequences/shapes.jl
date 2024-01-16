module Shapes
import ..Methods: start_time, end_time, add_TR

"""
    Shape(times, amplitudes; normalise=false)

Defines a [`Shape`](@ref) profile for an RF pulse or gradient profile.
The `Shape` is parametrised by the number of control points `N` and the type of the amplitude object `T`.
This amplitude type will be `Float64` for the amplitude and phase of the [`RFPulse`](@ref).
It will be `SVector{3, Float64}` for `MRGradients`.
"""
struct Shape{T}
    times :: Vector{Float64}
    amplitudes :: Vector{T}
    function Shape(times::AbstractVector{<:Number}, amplitudes::AbstractVector{T}) where {T}
        N = length(times)
        @assert N == length(amplitudes)
        if N == 0 || times[1] == times[end]
            return new{T}(Float64[], T[])
        end

        times = Float64.(times)

        # removing duplicate times
        prev_matched = false
        to_remove = Int[]
        for index in 1:N-1
            if times[index] == times[index + 1]
                if prev_matched
                    push!(to_remove, index)
                elseif index != 1
                    times[index] = prevfloat(times[index])
                end
                prev_matched = true
            elseif prev_matched
                times[index] = nextfloat(times[index])
                prev_matched = false
            end
        end
        to_keep = [i for i in 1:N if !(i in to_remove)]
        times = times[to_keep]
        amplitudes = amplitudes[to_keep]

        @assert all(times[2:end] .> times[1:end-1])
        new{T}(times, amplitudes)
    end
end

get_index(shape::Shape, time::Number) = time > end_time(shape) ? 0 : searchsortedlast(shape.times, time)
control_points(shape::Shape) = shape.times

Shape{T}() where {T} = Shape(Float64[], T[])

Base.length(shape::Shape) = length(shape.times)

"""
    sample(shape, time)

Returns the value of the shape at a given time using linear interpolation.
"""
function sample(shape::Shape{T}, time::Number) where {T}
    index = get_index(shape, time)
    if iszero(index)
        return zero(T)
    end
    dt0 = time - shape.times[index]
    if iszero(dt0)
        return shape.amplitudes[index]
    end
    dt1 = shape.times[index + 1] - time
    return (shape.amplitudes[index] * dt1 + shape.amplitudes[index + 1] * dt0) / (dt1 + dt0)
end

"""
    sample(shape, t1, t2)

Returns the average value of the shape between `t1` and `t2`. 
See [`sample_integral`](@ref) to get the full integral.
"""
sample(shape::Shape, t1::Number, t2::Number) = sample_integral(shape, t1, t2) / (t2 - t1)

"""
    sample_derivative(shape, time[, later_time])

Returns the derivative of the shape value at the specific time.
At control points the derivative is always assumed to be the slope of the next block.

If `later_time` is provided the average derivative between both timepoints is computed.
"""
function sample_derivative(shape::Shape{T}, time::Number) where {T}
    index = get_index(shape, time)
    if iszero(index) || index == length(shape)
        return zero(T)
    end
    return (shape.amplitudes[index + 1] - shape.amplitudes[index]) / (shape.times[index + 1] - shape.times[index])
end

sample_derivative(shape::Shape, t1::Number, t2::Number) = (sample(shape, t2) - sample(shape, t1)) / (t2 - t1)

"""
    sample_integral_step(shape, index, t0=`start of block`, t1=`end of block`)

Integrate the shape value from t0 to t1 assuming that both are within the single step (given by `index`).
This is a helper function for [`sample_integral`](@ref).
"""
function sample_integral_step(shape::Shape, index::Int, t0=nothing, t1=nothing)
    a0 = isnothing(t0) ? shape.amplitudes[index] : sample(shape, t0)
    a1 = isnothing(t1) ? shape.amplitudes[index+1] : sample(shape, t1)
    t0_real = isnothing(t0) ? shape.times[index] : t0
    t1_real = isnothing(t1) ? shape.times[index+1] : t1
    mean_amplitude = (a0 + a1) / 2
    return mean_amplitude * (t1_real - t0_real)
end

"""
    sample_integral(shape, t0, t1)

Integrate the shape value from t0 to t1.
"""
function sample_integral(shape::Shape{T}, t0::Number=-Inf, t1::Number=Inf) where {T}
    @assert t1 >= t0
    if (t0 > end_time(shape)) || (t1 < start_time(shape))
        return zero(T)
    end
    t0 = max(t0, start_time(shape))
    t1 = min(t1, end_time(shape))
    index0 = get_index(shape, t0)
    index1 = get_index(shape, t1)
    if index0 == index1
        return sample_integral_step(shape, index0, t0, t1)
    end
    total = sample_integral_step(shape, index0, t0, nothing)
    for index in index0+1:index1-1
        total += sample_integral_step(shape, index)
    end
    total += sample_integral_step(shape, index1, nothing, t1)
    return total
end

"""
    sample_integral(shape, t0, times)

Integrate the shape value from t0 to all the times in `times` (all assumed to be larger than `t0` and strictly increasing).
"""
function sample_integral(shape::Shape, t0::T, times::AbstractVector{T}) where {T <: Number}
    totals = zeros(Float64, length(times))
    current = zero(Float64)
    prev_time = t0
    for (index, t) in enumerate(times)
        current += sample_integral(shape, prev_time, t)
        totals[index] = current
    end
    return totals
end

add_TR(shape::Shape, TR::Number) = Shape(shape.times .+ TR, shape.amplitudes)

start_time(shape::Shape) = iszero(length(shape)) ? 0 : shape.times[1]
end_time(shape::Shape) = iszero(length(shape)) ? 0 : shape.times[end]


"""
    ShapePart(shape, t0, t1)

Represents a small part of a [`Shape`](@ref) between `t0` and `t1` during which the amplitude varies linearly.
This object is used to store the relevant part of the RF and gradient profiles within a single timestep (see [`SequencePart`](@ref)).

Like the full shape values can be extracted using [`sample`](@ref), [`sample_derivative`](@ref), [`sample_integral`](@ref).
The time passed on to these functions should be between 0 (for `t0`) and 1 (for `t1`).
"""
struct ShapePart{T}
    start :: T
    final :: T
    slope :: T
end

function ShapePart(shape::Shape{T}, t0::Number, t1::Number) where {T}
    tmean = (t0 + t1) / 2
    mean_value = sample(shape, tmean)
    slope = sample_derivative(shape, tmean)
    if any(isinf.(slope))
        slope = map(slope) do val
            isinf(val) ? prevfloat(val) : val
        end
    end
    return ShapePart{T}(
        (t0 - tmean) * slope .+ mean_value,
        (t1 - tmean) * slope .+ mean_value,
        slope * (t1 - t0),
    )
end

sample(part::ShapePart, time) = part.start + time * part.slope
sample(part::ShapePart, t1, t2) = sample(part, (t1 + t2) / 2)
sample_derivative(part::ShapePart, t1, t2=nothing) = part.slope
sample_integral(part::ShapePart, t1, t2) = sample(part, (t1 + t2) / 2) * (t2 - t1)
end