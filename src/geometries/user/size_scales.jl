module SizeScales
import ..Obstructions: ObstructionGroup, Walls, Cylinders, Annuli, Spheres, Mesh, isglobal
import ...Internal.Obstructions.Triangles: normal, curvature

"""
    size_scale(obstruction; ignore_user_value=false)

Determines the size scale of the given obstruction.

The user can set the size scale using `obstruction.size_scale=<value>`.
If set, this value will be returned (unless `ignore_user_value` is set to true).
"""
function size_scale(obstruction::ObstructionGroup; ignore_user_value=false)
    if ~ignore_user_value && ~isnothing(obstruction.size_scale.value)
        return obstruction.size_scale.value
    end
    return compute_size_scale(obstruction)
end


function compute_size_scale(obstruction::ObstructionGroup)
    repeats = isnothing(obstruction.repeats.value) ? Inf : minimum(obstruction.repeats.value)
    return min(
        repeats,
        compute_obstruction_size_scale(obstruction)
    )
end

function compute_obstruction_size_scale(w::Walls)
    if isglobal(w.position)
        p = [w.position.value]
    else
        p = sort(w.position.value)
    end
    if ~isnothing(w.repeats.value)
        push!(p, p[1] .+ w.repeats.value[1])
    end
    return minimum(p[2:end] - p[1:end-1]; init=Inf)
end
compute_obstruction_size_scale(c::Cylinders) = minimum(c.radius.value)
compute_obstruction_size_scale(a::Annuli) = min(
    minimum(a.inner.value),
    minimum(a.outer.value),
    minimum(a.outer.value .- a.inner.value),
)
compute_obstruction_size_scale(s::Spheres) = minimum(s.radius.value)
function compute_obstruction_size_scale(m::Mesh)
    c = m.n_obstructions > 1 ? curvature(m.triangles.value, m.vertices.value) : 0.
    return 1 / (2 * c)
end

end