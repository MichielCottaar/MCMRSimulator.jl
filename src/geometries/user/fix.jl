module Fix
import LinearAlgebra: I, transpose
import StaticArrays: SVector
import ...Internal: FixedObstructionGroup, FixedGeometry, Internal, get_quadrant, Grid, isinside_grid, FixedMesh
import ..Obstructions: ObstructionType, ObstructionGroup, Walls, Cylinders, Spheres, Annuli, Mesh, fields, isglobal
import ..SplitMesh: split_mesh

function fix(geometry::AbstractVector)
    result = FixedObstructionGroup[]
    for og in geometry
        append!(result, fix(og, length(result) + 1))
    end
    return Tuple(result)
end

function fix(geometry::ObstructionGroup, index::Int=1)
    for key in propertynames(geometry)
        field_value = getproperty(geometry, key)
        if field_value.field.required && any(isnothing.(field_value.value))
            error("Missing value for required field $field_value")
        end
    end
    result = fix_type(geometry, index)
    return result isa FixedObstructionGroup ? (result, ) : Tuple(result)
end


"""
    fix_type(group::ObstructionGroup, index)

Replace user [`ObstructionGroup`](@ref) objects with a [`FixedGeometry`](@ref).

This needs to be defined for each sub-type of [`ObstructionGroup`](@ref),
    e.g. [`Spheres`](@ref), [`Annuli`](@ref), [`Mesh`](@ref).

Shifts, collision properties, and volumetric and surface MRI properties can be added using [`apply_properties`](@ref).
"""
function fix_type end

function fix_type(walls::Walls, index::Int)
    base_obstructions = fill(Internal.Wall(), length(walls))
    apply_properties(walls, base_obstructions, index; surface="surface")
end

function fix_type(cylinders::Cylinders, index::Int)
    if isglobal(cylinders.radius)
        base_obstructions = fill(Internal.Cylinder(cylinders.radius.value), length(cylinders))
    else
        base_obstructions = [Internal.Cylinder(r) for r in cylinders.radius.value]
    end
    apply_properties(cylinders, base_obstructions, index; surface="surface", volume="volume")
end

function fix_type(spheres::Spheres, index::Int)
    if isglobal(spheres.radius)
        base_obstructions = fill(Internal.Sphere(spheres.radius.value), length(spheres))
    else
        base_obstructions = [Internal.Sphere(r) for r in spheres.radius.value]
    end
    apply_properties(spheres, base_obstructions, index; surface="surface", volume="volume")
end

function fix_type(annuli::Annuli, index::Int)
    if isglobal(annuli.inner)
        base_inner = fill(Internal.Cylinder(annuli.inner.value), length(annuli))
    else
        base_inner = [Internal.Cylinder(r) for r in annuli.inner.value]
    end
    inner_cylinders = apply_properties(annuli, base_inner, index; surface="inner_surface", volume="inner_volume")

    if isglobal(annuli.outer)
        base_outer = fill(Internal.Cylinder(annuli.outer.value), length(annuli))
    else
        base_outer = [Internal.Cylinder(r) for r in annuli.outer.value]
    end
    outer_cylinders = apply_properties(annuli, base_outer, index + 1; surface="outer_surface", volume="outer_volume")
    return (inner_cylinders, outer_cylinders)
end

function fix_type(mesh::Mesh, index::Int)
    [fix_type_single(m, index + i - 1) for (i, m) in enumerate(split_mesh(mesh))]
end

function fix_type_single(mesh::Mesh, index::Int)
    base_obstructions = [Internal.IndexTriangle(index) for index in mesh.triangles.value]
    if ~mesh.save_memory.value
        base_obstructions = [Internal.FullTriangle(it, SVector{3}.(mesh.vertices.value)) for it in base_obstructions]
    end
    apply_properties(mesh, base_obstructions, index; surface="surface", volume="volume")
end

"""
    apply_properties(user_obstructions, internal_obstructions, index; surface, volume)

Apply the generic obstruction properties defined in `user_obstruction` to `internal_obstructions`.
This function applies (if appropriate):
- positional shifts
- volumetic and surface-bound MRI properties
- collision properties
- rotation
- repeats
- mesh vertices
"""
function apply_properties(user_obstructions::ObstructionGroup, internal_obstructions::Vector{<:Internal.FixedObstruction}, index::Int; surface=nothing, volume=nothing)
    # apply shifts
    if hasproperty(user_obstructions, :position)
        shifts = isglobal(user_obstructions.position) ? fill(user_obstructions.position.value, length(user_obstructions)) : user_obstructions.position.value
        if hasproperty(user_obstructions, :repeats) && ~isnothing(user_obstructions.repeats.value)
            shifts = [mod.(s .+ user_obstructions.repeats.value/2, user_obstructions.repeats.value) .- user_obstructions.repeats.value/2 for s in shifts]
        end
        internal_obstructions = Internal.Shift.(internal_obstructions, shifts)
    end

    fix_array_type(other) = other
    function fix_array_type(arr::Vector{Union{Nothing, T}}) where {T}
        if ~any(isnothing.(arr))
            return Vector{T}(arr)
        else
            return arr
        end
    end

    # get volumetric MRI properties
    if ~isnothing(volume)
        get_volume(s) = fix_array_type(getproperty(user_obstructions, Symbol(String(s) * "_" * String(volume))).value)
        symbols = filter(s->~isnothing(get_volume(s)), [:R1, :R2, :off_resonance])
        values = [get_volume(s) for s in symbols]
        volume = NamedTuple{Tuple(symbols)}(Tuple(values))
    else
        volume = NamedTuple{()}(())
    end

    # get surface properties (bound MRI & collision)
    if ~isnothing(surface)
        get_surface(s) = fix_array_type(getproperty(user_obstructions, Symbol(String(s) * "_" * String(surface))).value)
        symbols = filter(s->~isnothing(get_surface(s)), [:R1, :R2, :off_resonance, :permeability, :density, :dwell_time, :relaxivity])
        values = [get_surface(s) for s in symbols]
        updated_symbols = replace(symbols, :density=>:surface_density, :relaxivity=>:surface_relaxivity)
        surface = NamedTuple{Tuple(updated_symbols)}(Tuple(values))
    else
        surface = NamedTuple{()}(())
    end

    rotation = user_obstructions.rotation.value

    if eltype(internal_obstructions) <: Internal.IndexTriangle
        vertices = SVector{3, Float64}.(user_obstructions.vertices.value)
    else
        vertices = SVector{3, Float64}[]
    end
    bounding_boxes = map(o->Internal.BoundingBox(o, vertices), internal_obstructions)
    grid = Grid(bounding_boxes, user_obstructions.grid_resolution.value, user_obstructions.repeats.value)
    
    result = FixedObstructionGroup(
        internal_obstructions,
        index,
        rotation,
        grid,
        bounding_boxes,
        volume,
        surface,
        vertices
    )
    if result isa FixedMesh
        grid = Grid(bounding_boxes, user_obstructions.grid_resolution.value, user_obstructions.repeats.value; isinside=isinside_grid(result))
        result = FixedObstructionGroup(
            internal_obstructions,
            index,
            rotation,
            grid,
            bounding_boxes,
            volume,
            surface,
            vertices
        )
    end
    return result
end


end