"""Thread-safe metric registry for opt-in diagnostics."""
module DebugCounters

export AbstractMetric, FrequencyTable, SummaryMeasures, Histogram, register!, record!, reset!, snapshot, report

abstract type AbstractMetric end

"""Exact occurrence counts for each observed integer value."""
mutable struct FrequencyTable{T<:Integer} <: AbstractMetric
    counts::Dict{T, Int}
end

FrequencyTable(::Type{T}) where {T<:Integer} = FrequencyTable{T}(Dict{T, Int}())

"""Online count, extrema, mean, and sample variance for a registered metric."""
mutable struct SummaryMeasures <: AbstractMetric
    count::Int
    minimum::Float64
    maximum::Float64
    mean::Float64
    m2::Float64
end

SummaryMeasures() = SummaryMeasures(0, Inf, -Inf, 0.0, 0.0)

"""Fixed-edge histogram, with separate counts for values outside its edges."""
mutable struct Histogram <: AbstractMetric
    edges::Vector{Float64}
    counts::Vector{Int}
    underflow::Int
    overflow::Int
end

function Histogram(edges)
    converted = Float64.(collect(edges))
    length(converted) >= 2 || throw(ArgumentError("histogram bins must contain at least two edges"))
    all(isfinite, converted) || throw(ArgumentError("histogram edges must be finite"))
    all(diff(converted) .> 0) || throw(ArgumentError("histogram edges must be strictly increasing"))
    Histogram(converted, zeros(Int, length(converted) - 1), 0, 0)
end

const _lock = ReentrantLock()
const _metrics = Dict{Symbol, Tuple{DataType, AbstractMetric}}()

_metric_configuration(::Type{T}; histogram=false, bins=nothing) where {T<:Real} =
    _metric_configuration(T, histogram, bins)

function _metric_configuration(::Type{T}, histogram::Bool, bins) where {T<:Real}
    if histogram
        isnothing(bins) && throw(ArgumentError("histogram=true requires explicit bin edges"))
        return T, Histogram(bins)
    end
    isnothing(bins) || throw(ArgumentError("bins can only be provided with histogram=true"))
    T <: Integer && return T, FrequencyTable(T)
    T, SummaryMeasures()
end

"""
    register!(name::Symbol, T::Type{<:Real}; histogram=false, bins=nothing)

Register a real-valued metric. Integer types use exact frequency tables by
default; other real types collect summary measures (count, minimum, maximum,
mean, and sample variance). Set `histogram=true` and provide explicit `bins`
edges to collect a histogram with underflow and overflow counts.
"""
function register!(name::Symbol, ::Type{T}; histogram=false, bins=nothing) where {T<:Real}
    isconcretetype(T) || throw(ArgumentError("registered metric type must be concrete"))
    histogram isa Bool || throw(ArgumentError("histogram must be a Bool"))
    value_type, metric = _metric_configuration(T; histogram, bins)
    lock(_lock) do
        if haskey(_metrics, name)
            existing_type, existing_metric = _metrics[name]
            compatible = existing_type === value_type &&
                typeof(existing_metric) === typeof(metric) &&
                (!(metric isa Histogram) || existing_metric.edges == metric.edges)
            compatible || throw(ArgumentError("metric :$name is already registered with a different configuration"))
            return existing_metric
        end
        _metrics[name] = (value_type, metric)
        metric
    end
end

function _convert_value(::Type{T}, value::Real) where {T<:Real}
    converted = try
        convert(T, value)
    catch
        throw(ArgumentError("value $value cannot be converted to registered metric type $T"))
    end
    isequal(converted, value) || throw(ArgumentError(
        "value $value cannot be represented exactly by registered metric type $T",
    ))
    isfinite(converted) || throw(ArgumentError("metric values must be finite"))
    converted
end

function _update!(metric::FrequencyTable{T}, value::T) where {T<:Integer}
    metric.counts[value] = get(metric.counts, value, 0) + 1
    metric
end

function _update!(metric::SummaryMeasures, real_value::Real)
    value = Float64(real_value)
    metric.count += 1
    metric.minimum = min(metric.minimum, value)
    metric.maximum = max(metric.maximum, value)
    delta = value - metric.mean
    metric.mean += delta / metric.count
    metric.m2 += delta * (value - metric.mean)
    metric
end

function _update!(metric::Histogram, real_value::Real)
    value = Float64(real_value)
    if value < first(metric.edges)
        metric.underflow += 1
    elseif value > last(metric.edges)
        metric.overflow += 1
    else
        bin = value == last(metric.edges) ? length(metric.counts) : searchsortedlast(metric.edges, value)
        metric.counts[bin] += 1
    end
    metric
end

"""Record a value in a previously registered metric."""
function record!(name::Symbol, value::Real)
    lock(_lock) do
        haskey(_metrics, name) || throw(ArgumentError("metric :$name has not been registered"))
        value_type, metric = _metrics[name]
        _update!(metric, _convert_value(value_type, value))
    end
    nothing
end

"""Remove all metric registrations and recorded values."""
function reset!()
    lock(_lock) do
        empty!(_metrics)
    end
    nothing
end

_snapshot(metric::FrequencyTable) = copy(metric.counts)

function _snapshot(metric::SummaryMeasures)
    (
        count=metric.count,
        minimum=metric.count == 0 ? nothing : metric.minimum,
        maximum=metric.count == 0 ? nothing : metric.maximum,
        mean=metric.count == 0 ? nothing : metric.mean,
        variance=metric.count < 2 ? nothing : metric.m2 / (metric.count - 1),
    )
end

_snapshot(metric::Histogram) = (
    edges=copy(metric.edges),
    counts=copy(metric.counts),
    underflow=metric.underflow,
    overflow=metric.overflow,
)

"""Return an independent snapshot of all registered metrics."""
function snapshot()
    lock(_lock) do
        Dict(name => _snapshot(metric) for (name, (_, metric)) in _metrics)
    end
end

"""Print summaries of all registered metrics."""
function report(io::IO=stdout)
    values = snapshot()
    for name in sort!(collect(keys(values)); by=string)
        println(io, name, ": ", values[name])
    end
    nothing
end

end
