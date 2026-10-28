"""Metric collectors and thread-safe tables for opt-in diagnostics."""
module DebugCounters

export AbstractMetric, DebugTable, FrequencyCounter, FrequencyTable,
    SummaryMeasures, Histogram, record!, reset!, snapshot, report

abstract type AbstractMetric end

"""Exact occurrence counts for integer values."""
mutable struct FrequencyCounter{T<:Integer} <: AbstractMetric
    counts::Dict{T, Int}
end

FrequencyCounter(::Type{T}) where {T<:Integer} = FrequencyCounter{T}(Dict{T, Int}())

# Compatibility name used by the current Grid debug state.
const FrequencyTable = FrequencyCounter

function record!(metric::FrequencyCounter{T}, value::Integer) where {T<:Integer}
    converted = try
        convert(T, value)
    catch
        throw(ArgumentError("value $value cannot be represented by frequency-counter key type $T"))
    end
    converted == value || throw(ArgumentError(
        "value $value cannot be represented exactly by frequency-counter key type $T",
    ))
    metric.counts[converted] = get(metric.counts, converted, 0) + 1
    nothing
end

record!(::FrequencyCounter, value::Real) = throw(ArgumentError(
    "frequency counters accept only integer values, got $(typeof(value))",
))

reset!(metric::FrequencyCounter) = (empty!(metric.counts); nothing)

function report(io::IO, metric::FrequencyCounter)
    show(io, metric.counts)
    nothing
end

"""Online count, extrema, mean, and sample variance."""
mutable struct SummaryMeasures <: AbstractMetric
    count::Int
    minimum::Float64
    maximum::Float64
    mean::Float64
    m2::Float64
end

SummaryMeasures() = SummaryMeasures(0, Inf, -Inf, 0.0, 0.0)

function record!(metric::SummaryMeasures, value::Real)
    value = Float64(value)
    isfinite(value) || throw(ArgumentError("metric values must be finite"))
    metric.count += 1
    metric.minimum = min(metric.minimum, value)
    metric.maximum = max(metric.maximum, value)
    delta = value - metric.mean
    metric.mean += delta / metric.count
    metric.m2 += delta * (value - metric.mean)
    nothing
end

function reset!(metric::SummaryMeasures)
    metric.count = 0
    metric.minimum = Inf
    metric.maximum = -Inf
    metric.mean = 0.0
    metric.m2 = 0.0
    nothing
end

function report(io::IO, metric::SummaryMeasures)
    show(io, (
        count=metric.count,
        minimum=metric.count == 0 ? nothing : metric.minimum,
        maximum=metric.count == 0 ? nothing : metric.maximum,
        mean=metric.count == 0 ? nothing : metric.mean,
        variance=metric.count < 2 ? nothing : metric.m2 / (metric.count - 1),
    ))
    nothing
end

"""Fixed-edge histogram with separate underflow and overflow counts."""
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

function record!(metric::Histogram, real_value::Real)
    value = Float64(real_value)
    isfinite(value) || throw(ArgumentError("metric values must be finite"))
    if value < first(metric.edges)
        metric.underflow += 1
    elseif value > last(metric.edges)
        metric.overflow += 1
    else
        bin = value == last(metric.edges) ? length(metric.counts) : searchsortedlast(metric.edges, value)
        metric.counts[bin] += 1
    end
    nothing
end

reset!(metric::Histogram) = begin
    fill!(metric.counts, 0)
    metric.underflow = 0
    metric.overflow = 0
    nothing
end

function report(io::IO, metric::Histogram)
    show(io, (
        edges=metric.edges,
        counts=metric.counts,
        underflow=metric.underflow,
        overflow=metric.overflow,
    ))
    nothing
end

"""Named, described collectors protected by one table-level lock."""
struct DebugTable
    order::Vector{Symbol}
    metrics::Dict{Symbol, Tuple{AbstractMetric, String}}
    lock::ReentrantLock
end

"""
    DebugTable(metrics)

Create a table from `(name, description, metric)` entries. Each name must be
unique. Access its collectors through `record!`, `reset!`, `snapshot`, and
`report` on the table so its single lock protects every operation.
"""
function DebugTable(metrics)
    order = Symbol[]
    metric_dict = Dict{Symbol, Tuple{AbstractMetric, String}}()
    for (name, description, metric) in metrics
        name isa Symbol || throw(ArgumentError("metric names must be Symbols"))
        description isa AbstractString || throw(ArgumentError("metric descriptions must be strings"))
        metric isa AbstractMetric || throw(ArgumentError("metric $name must be an AbstractMetric"))
        haskey(metric_dict, name) && throw(ArgumentError("duplicate metric name :$name"))
        push!(order, name)
        metric_dict[name] = (metric, String(description))
    end
    DebugTable(order, metric_dict, ReentrantLock())
end

"""Record a value in the named collector under the table lock."""
function record!(table::DebugTable, name::Symbol, value)
    lock(table.lock) do
        haskey(table.metrics, name) || throw(KeyError(name))
        record!(table.metrics[name][1], value)
    end
    nothing
end

"""Reset every collector in the table under the table lock."""
function reset!(table::DebugTable)
    lock(table.lock) do
        for name in table.order
            reset!(table.metrics[name][1])
        end
    end
    nothing
end

"""Return an independent snapshot of one named metric."""
function snapshot(table::DebugTable, name::Symbol)
    lock(table.lock) do
        haskey(table.metrics, name) || throw(KeyError(name))
        deepcopy(table.metrics[name][1])
    end
end

"""Return independent snapshots of all metrics, keyed by metric name."""
function snapshot(table::DebugTable)
    lock(table.lock) do
        Dict(name => deepcopy(table.metrics[name][1]) for name in table.order)
    end
end

"""Print a table's descriptions and current metric snapshots."""
function report(io::IO, table::DebugTable)
    metrics = snapshot(table)
    for name in table.order
        println(io, "# ", name, ": ", table.metrics[name][2])
        report(io, metrics[name])
        println(io)
    end
    nothing
end

report(table::DebugTable; io::IO=stdout) = report(io, table)

end
