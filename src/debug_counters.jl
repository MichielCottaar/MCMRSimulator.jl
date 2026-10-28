"""Thread-safe metric registry for opt-in diagnostics."""
module DebugCounters

export AbstractMetric, FrequencyTable, SummaryMeasures, Histogram, register!, record!, reset!, snapshot, report

abstract type AbstractMetric end

"""Exact occurrence counts for each observed integer value."""
mutable struct FrequencyTable{T<:Integer} <: AbstractMetric
    counts::Dict{T, Int}
    lock::ReentrantLock
end

FrequencyTable(::Type{T}) where {T<:Integer} =
    FrequencyTable{T}(Dict{T, Int}(), ReentrantLock())

"""Record one occurrence in a local frequency table."""
function record!(table::FrequencyTable{T}, value::Integer) where {T<:Integer}
    converted = try
        convert(T, value)
    catch
        throw(ArgumentError("value $value cannot be represented by frequency-table key type $T"))
    end
    converted == value || throw(ArgumentError(
        "value $value cannot be represented exactly by frequency-table key type $T",
    ))
    lock(table.lock) do
        table.counts[converted] = get(table.counts, converted, 0) + 1
    end
    nothing
end

"""Return a copy of a local frequency table's counts."""
snapshot(table::FrequencyTable) = lock(table.lock) do
    copy(table.counts)
end

function report(io::IO, table::FrequencyTable)
    show(io, snapshot(table))
    nothing
end
report(table::FrequencyTable; io::IO=stdout) = report(io, table)

"""Online count, extrema, mean, and sample variance for a registered metric."""
mutable struct SummaryMeasures <: AbstractMetric
    count::Int
    minimum::Float64
    maximum::Float64
    mean::Float64
    m2::Float64
    lock::ReentrantLock
end

SummaryMeasures() = SummaryMeasures(0, Inf, -Inf, 0.0, 0.0, ReentrantLock())

"""Fixed-edge histogram, with separate counts for values outside its edges."""
mutable struct Histogram <: AbstractMetric
    edges::Vector{Float64}
    counts::Vector{Int}
    underflow::Int
    overflow::Int
    lock::ReentrantLock
end

function Histogram(edges)
    converted = Float64.(collect(edges))
    length(converted) >= 2 || throw(ArgumentError("histogram bins must contain at least two edges"))
    all(isfinite, converted) || throw(ArgumentError("histogram edges must be finite"))
    all(diff(converted) .> 0) || throw(ArgumentError("histogram edges must be strictly increasing"))
    Histogram(converted, zeros(Int, length(converted) - 1), 0, 0, ReentrantLock())
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
    lock(metric.lock) do
        metric.counts[value] = get(metric.counts, value, 0) + 1
    end
    metric
end

function _update!(metric::SummaryMeasures, real_value::Real)
    value = Float64(real_value)
    isfinite(value) || throw(ArgumentError("metric values must be finite"))
    lock(metric.lock) do
        metric.count += 1
        metric.minimum = min(metric.minimum, value)
        metric.maximum = max(metric.maximum, value)
        delta = value - metric.mean
        metric.mean += delta / metric.count
        metric.m2 += delta * (value - metric.mean)
    end
    metric
end

function _update!(metric::Histogram, real_value::Real)
    value = Float64(real_value)
    isfinite(value) || throw(ArgumentError("metric values must be finite"))
    lock(metric.lock) do
        if value < first(metric.edges)
            metric.underflow += 1
        elseif value > last(metric.edges)
            metric.overflow += 1
        else
            bin = value == last(metric.edges) ? length(metric.counts) : searchsortedlast(metric.edges, value)
            metric.counts[bin] += 1
        end
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

"""Record a value directly in a local summary or histogram collector."""
record!(metric::Union{SummaryMeasures, Histogram}, value::Real) = _update!(metric, value)

"""Remove all metric registrations and recorded values."""
function reset!()
    lock(_lock) do
        empty!(_metrics)
    end
    nothing
end

function snapshot(metric::SummaryMeasures)
    lock(metric.lock) do
        (
            count=metric.count,
            minimum=metric.count == 0 ? nothing : metric.minimum,
            maximum=metric.count == 0 ? nothing : metric.maximum,
            mean=metric.count == 0 ? nothing : metric.mean,
            variance=metric.count < 2 ? nothing : metric.m2 / (metric.count - 1),
        )
    end
end

function snapshot(metric::Histogram)
    lock(metric.lock) do
        (
            edges=copy(metric.edges),
            counts=copy(metric.counts),
            underflow=metric.underflow,
            overflow=metric.overflow,
        )
    end
end

"""Return an independent snapshot of all registered metrics."""
function snapshot()
    lock(_lock) do
        Dict(name => snapshot(metric) for (name, (_, metric)) in _metrics)
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
