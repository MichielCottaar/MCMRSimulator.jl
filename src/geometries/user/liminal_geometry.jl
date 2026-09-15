module LiminalGeometries

export LiminalGeometry

"""
    LiminalGeometry(geometries; extracellular_fraction)

Create a statistical collection of cell geometry templates. `geometries` should
be a vector of `(number_fraction, geometry)` tuples. The number fractions give
the relative abundance of each cell type and are normalized internally;
`extracellular_fraction` gives the extracellular volume fraction.

See the [Liminal geometry](@ref liminal_geometry) chapter for an overview and
usage examples.
"""
struct LiminalGeometry
    geometries::Vector{Tuple{Float64, Any}}
    extracellular_fraction::Float64
end

function LiminalGeometry(
    geometries::AbstractVector;
    extracellular_fraction::Real,
)
    isempty(geometries) && throw(ArgumentError("at least one cell geometry is required"))
    isfinite(extracellular_fraction) && 0 <= extracellular_fraction <= 1 ||
        throw(ArgumentError("extracellular_fraction must be between 0 and 1"))

    entries = Tuple{Float64, Any}[]
    for entry in geometries
        entry isa Tuple && length(entry) == 2 ||
            throw(ArgumentError("each cell geometry must be a (number_fraction, geometry) tuple"))
        fraction, geometry = entry
        fraction isa Real && isfinite(fraction) && fraction > 0 ||
            throw(ArgumentError("cell number fractions must be finite and positive"))
        push!(entries, (Float64(fraction), geometry))
    end
    LiminalGeometry(entries, Float64(extracellular_fraction))
end

LiminalGeometry(geometries::Tuple; kwargs...) =
    LiminalGeometry(collect(geometries); kwargs...)

LiminalGeometry(; geometries, kwargs...) =
    LiminalGeometry(geometries; kwargs...)

end
