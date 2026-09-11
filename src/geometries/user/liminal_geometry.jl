module LiminalGeometries

export LiminalGeometry

"""User-facing collection of cell geometry templates and volume fractions."""
struct LiminalGeometry
    geometries::Vector{Tuple{Float64, Any}}
    extracellular_fraction::Float64
end

function LiminalGeometry(
    geometries::AbstractVector;
    extracellular_fraction::Real=0.,
)
    isempty(geometries) && throw(ArgumentError("at least one cell geometry is required"))
    isfinite(extracellular_fraction) && 0 <= extracellular_fraction <= 1 ||
        throw(ArgumentError("extracellular_fraction must be between 0 and 1"))

    entries = Tuple{Float64, Any}[]
    for entry in geometries
        entry isa Tuple && length(entry) == 2 ||
            throw(ArgumentError("each cell geometry must be a (volume_fraction, geometry) tuple"))
        fraction, geometry = entry
        fraction isa Real && isfinite(fraction) && fraction > 0 ||
            throw(ArgumentError("cell volume fractions must be finite and positive"))
        push!(entries, (Float64(fraction), geometry))
    end
    LiminalGeometry(entries, Float64(extracellular_fraction))
end

LiminalGeometry(geometries::Tuple; kwargs...) =
    LiminalGeometry(collect(geometries); kwargs...)

LiminalGeometry(; geometries, kwargs...) =
    LiminalGeometry(geometries; kwargs...)

end
