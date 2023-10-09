module FixSusceptibility
import StaticArrays: SVector
import LinearAlgebra: transpose
import ...Internal.Susceptibility: FixedSusceptibility, ParentSusceptibility, BaseSusceptibility, CylinderSusceptibility, AnnulusSusceptibility, TriangleSusceptibility
import ...Internal: Grid, BoundingBox, FullTriangle, radius
import ..Obstructions: ObstructionType, ObstructionGroup, Walls, Cylinders, Spheres, Annuli, Mesh, fields, isglobal

"""
    fix_susceptibility(geometry)

Create a new [`FixedSusceptibility`](@ref) based on the user-defined geometry settings.
"""
function fix_susceptibility(geometry::AbstractVector)
    result = ParentSusceptibility[]
    for g in geometry
        addition = fix_susceptibility_type(g)
        if ~isnothing(addition)
            push!(result, addition)
        end
    end
    return Tuple(result)
end

fix_susceptibility(group::ObstructionGroup) = fix_susceptibility([group])


fix_susceptibility_type(unknown_group::ObstructionGroup) = nothing

function fix_susceptibility_type(group::Cylinders)
    if all(isone.(group.g_ratio.value))
        return nothing
    end
    b0_field = group.rotation.value[3, :]
    f(a...) = CylinderSusceptibility(a..., b0_field)
    base = f.(group.radius.value, group.g_ratio.value, group.susceptibility_iso.value, group.susceptibility_aniso.value)
    if ~(base isa Vector)
        base = fill(base, length(group))
    end
    add_parent(group, base)
end

function fix_susceptibility_type(group::Annuli)
    if ~any(group.myelin.value)
        return nothing
    end
    b0_field = group.rotation.value[3, :]
    f(a...) = AnnulusSusceptibility(a..., b0_field)
    base = f.(group.inner.value, group.outer.value, group.susceptibility_iso.value, group.susceptibility_aniso.value)
    if ~(base isa Vector)
        base = fill(base, length(group))
    end
    add_parent(group, base; radius_symbol=:inner)
end

function fix_susceptibility_type(group::Mesh)
    if ~any(group.myelin.value)
        return nothing
    end
    b0_field = group.rotation.value[3, :]
    res = TriangleSusceptibility[]
    radii = Float64[]
    positions = SVector{3, Float64}[]
    for (index, (i1, i2, i3)) in enumerate(group.triangles.value)
        ft = FullTriangle(group.vertices.value[i1], group.vertices.value[i2], group.vertices.value[i3])
        push!(radii, radius(ft))
        push!(positions, ft.a)
        push!(res, TriangleSusceptibility(ft, group[index].susceptibility_iso, group[index].susceptibility_aniso, b0_field))
    end
    return add_parent(group, res; positions=positions, radii=radii)
end

function add_parent(user::ObstructionGroup, internal::AbstractVector{<:BaseSusceptibility{N}}; positions=nothing, radii=nothing, radius_symbol=:radius) where {N}
    if isnothing(positions)
        positions = isglobal(user.position) ? fill(SVector{N}(user.position.value), length(internal)) : SVector{N}.(user.position.value)
    end

    if isnothing(radii)
        user_radius = getproperty(user, radius_symbol)
        radii = isglobal(user_radius) ? fill(user_radius.value, length(internal)) : user_radius.value
    end

    repeats = isnothing(user.repeats.value) ? nothing : SVector{N}(user.repeats.value)
    bbs = map((p, r) -> BoundingBox(p .- (r + user.lorentz_radius.value), p .+ (r + user.lorentz_radius.value)), positions, radii)
    grid = Grid(
        bbs,
        user.grid_resolution.value,
        repeats
    )
    ParentSusceptibility{length(internal), N, eltype(internal), typeof(repeats), N * 3}(
        internal,
        grid,
        Tuple{SVector{N, Float64}, Float64}[(p, r) for (p, r) in zip(positions, radii)],
        transpose(user.rotation.value),
        isnothing(repeats) ? nothing : repeats ./ 2,
        user.lorentz_radius.value,
        maximum(radii)
    )
end

end