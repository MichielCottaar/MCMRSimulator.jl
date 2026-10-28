"""
    DebugCounters

Thread-safe, per-instance metrics for opt-in diagnostics.

Call [`register!`](@ref) once for each instrumented instance type, providing
named metric templates and descriptions. Metrics are instances of an
[`AbstractMetric`](@ref) subtype:

- [`FrequencyCounter`](@ref) records exact integer-value frequencies and is
  constructed with its key type, for example `FrequencyCounter(Int)`.
- [`Counter`](@ref) accumulates signed integer increments.
- [`SummaryMeasures`](@ref) tracks count, extrema, mean, and sample variance.
- [`Histogram`](@ref) records values into explicitly specified edges, with
  separate underflow and overflow counts.

For example:

```julia
register!(MyInstrumentedType, [
    (:cells_visited, "Cells visited per query", FrequencyCounter(Int)),
    (:step_length, "Step length", SummaryMeasures()),
    (:step_distribution, "Step-length distribution", Histogram(0:1:10)),
])

record!(instance, :cells_visited, 3)
report(instance)
```

`register!` stores metric templates for the supplied type. On the first
`record!`, `snapshot`, or `report` for an instance, a private [`DebugTable`](@ref)
is created with copies of those templates. If registrations exist for both a
type and its supertype, the most-specific registration is used.

`record!(instance, name, value)`, `reset!(instance)`, `snapshot(instance)`, and
`report(instance)` operate on that instance's table under its lock. `reset!()`
without an argument discards all per-instance tables while preserving type
registrations. Snapshots contain independent copies of the metrics.
"""
module DebugCounters

export AbstractMetric, DebugTable, FrequencyCounter, FrequencyTable,
    Counter, SummaryMeasures, Histogram, register!, record!, reset!, snapshot, report

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

"""An integer accumulator updated by signed integer increments."""
mutable struct Counter{T<:Integer} <: AbstractMetric
    count::T
end

Counter(::Type{T}=Int) where {T<:Integer} = Counter{T}(zero(T))

function record!(metric::Counter{T}, increment::Integer) where {T<:Integer}
    converted = try
        convert(T, increment)
    catch
        throw(ArgumentError("increment $increment cannot be represented by counter type $T"))
    end
    converted == increment || throw(ArgumentError(
        "increment $increment cannot be represented exactly by counter type $T",
    ))
    metric.count = Base.checked_add(metric.count, converted)
    nothing
end

record!(::Counter, increment::Real) = throw(ArgumentError(
    "counters accept only integer increments, got $(typeof(increment))",
))

reset!(metric::Counter{T}) where {T<:Integer} = (metric.count = zero(T); nothing)

function report(io::IO, metric::Counter)
    show(io, metric.count)
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

const _registry_lock = ReentrantLock()
const _registered_tables = Dict{Any, DebugTable}()
const _instance_tables = Dict{UInt, DebugTable}()

function _metric_definitions(table::DebugTable)
    lock(table.lock) do
        [
            (name, table.metrics[name][2], deepcopy(table.metrics[name][1]))
            for name in table.order
        ]
    end
end

"""
    register!(T::Type, metrics)

Register `(name, description, metric)` definitions for instances of `T`.
Each instance receives its own lazily-created `DebugTable` with fresh copies of
these metric templates. A more-specific registered type takes precedence over
its registered supertypes.
"""
function register!(::Type{T}, metrics) where {T}
    table_template = DebugTable(metrics)
    definitions = _metric_definitions(table_template)
    lock(_registry_lock) do
        haskey(_registered_tables, T) && throw(ArgumentError(
            "debug metrics are already registered for $T",
        ))
        _registered_tables[T] = DebugTable(definitions)
    end
    nothing
end

function _registered_table_type(instance)
    instance_type = typeof(instance)
    matches = [registered for registered in keys(_registered_tables) if instance_type <: registered]
    isempty(matches) && throw(ArgumentError(
        "no debug metrics registered for instance type $instance_type",
    ))
    most_specific = [
        registered for registered in matches
        if !any(other -> other !== registered && other <: registered, matches)
    ]
    length(most_specific) == 1 || throw(ArgumentError(
        "ambiguous debug metric registrations for instance type $instance_type: $most_specific",
    ))
    only(most_specific)
end

function _instance_table(instance)
    instance_id = UInt(objectid(instance))
    lock(_registry_lock) do
        get!(_instance_tables, instance_id) do
            registered_type = _registered_table_type(instance)
            DebugTable(_metric_definitions(_registered_tables[registered_type]))
        end
    end
end

"""Record a metric value for this instance, creating its table on first use."""
record!(instance, name::Symbol, value) = record!(_instance_table(instance), name, value)

"""Reset the metrics associated with this instance."""
reset!(instance) = reset!(_instance_table(instance))

"""Return independent metric copies associated with this instance."""
snapshot(instance, name::Symbol) = snapshot(_instance_table(instance), name)
snapshot(instance) = snapshot(_instance_table(instance))

"""Report the metrics associated with this instance."""
report(io::IO, instance) = report(io, _instance_table(instance))
report(instance; io::IO=stdout) = report(io, instance)

"""Discard all per-instance tables while retaining registered metric definitions."""
function reset!()
    lock(_registry_lock) do
        empty!(_instance_tables)
    end
    nothing
end

end
