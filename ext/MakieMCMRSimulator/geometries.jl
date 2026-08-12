module Geometries
using Makie
import StaticArrays: SVector
import LinearAlgebra: cross, ⋅, norm
import Colors
import GeometryBasics
import MCMRSimulator.Plot: PlotPlane, GeometryLike, project_geometry, plot_geometry, plot_geometry!
import MCMRSimulator.Geometries.Internal: FixedGeometry, geometry_mesh
import MCMRSimulator.Geometries: ObstructionGroup, fix, Mesh, Cylinders, Spheres, Walls

@recipe Plot_Geometry (plot_plane::Union{Nothing, PlotPlane}, geometry::GeometryLike) begin
    "Set the color of the lines (2D) or patches (3D). In 2D it is set to the :black by default. In 3D each individual obstruction is by default plotted in a different, distinguishable color."
    color=Makie.automatic
    "Set the transparancy in a 3D plot (0 being fully transparent and 1 fully opague)."
    alpha=1.
    "Set the linewidth in 2D plots."
    linewidth=@inherit linewidth
    "Set the linestyle in 2D plots."
    linestyle=nothing
    "Set the line colour in 2D plots"
    linecolor=@inherit linecolor
    "Number of samples in mesh used to plot cylinders (default: 100) and spheres (default: 1000) in 3D plot."
    nsamples=Makie.automatic
    "Size to plot in μm of infinite walls and cylinders in 3D plot."
    height=1.


    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_shading_attributes()...
end

Makie.argument_names(::Type{<: Plot_Geometry}, N) = (:plot_plane, :geometry)

function Makie.plot!(scene::Plot_Geometry{<:Tuple{<:PlotPlane, <:GeometryLike}})
    map!(scene.attributes, [:color], :final_color) do color
        return color == Makie.automatic ? :black : color
    end
    map!(scene.attributes, [:plot_plane, :geometry], :line) do plot_plane, geometry
        return project_geometry(plot_plane, geometry)
    end

    Makie.lines!(scene, scene.attributes, scene.line; color=scene.final_color)
    scene
end

Makie.plottype(::PlotPlane, ::GeometryLike) = Plot_Geometry

Makie.convert_arguments(::Type{Plot_Geometry}, geometry::GeometryLike) = (nothing, geometry)


fix_and_mesh(geom::FixedGeometry; height, nsamples) = geom
fix_and_mesh(geom::Mesh; height, nsamples) = fix(geom)
fix_and_mesh(geom::Cylinders; height, nsamples) = fix(Mesh(geom; height=height, nsamples=nsamples == Makie.automatic ? 100 : nsamples))
fix_and_mesh(geom::Spheres; height, nsamples) = fix(Mesh(geom; nsamples=nsamples == Makie.automatic ? 1000 : nsamples))
fix_and_mesh(geom::Walls; height, nsamples) = fix(Mesh(geom; height=height))
fix_and_mesh(geom::ObstructionGroup; height, nsamples) = fix(Mesh(geom))

function Makie.plot!(scene::Plot_Geometry{<:Tuple{<:Nothing, <:GeometryLike}})
    Makie.register_computation!(scene.attributes, [:geometry, :color, :height, :nsamples], [:vertices, :triangles, :mesh_color]) do inputs, changed, cached
        mesh = fix_and_mesh(inputs[:geometry]; height=inputs[:height], nsamples=inputs[:nsamples])
        mesh_data = geometry_mesh(mesh)
        mesh_color_arr = Colors.distinguishable_colors(length(mesh_data))
        vertices = GeometryBasics.Point{3, Float64}[]
        triangles = GeometryBasics.TriangleFace{Int}[]
        colors = Any[]

        for (index, group) in enumerate(mesh_data)
            append!(triangles, [GeometryBasics.TriangleFace{Int}(t .+ length(vertices)) for t in group.triangles])
            append!(vertices, GeometryBasics.Point{3, Float64}.(group.vertices))
            patch_color = inputs[:color] == Makie.automatic ? mesh_color_arr[index] : inputs[:color]
            append!(colors, [patch_color for _ in group.vertices])
        end
        return (vertices, triangles, colors)
    end

    Makie.mesh!(scene, scene.attributes, scene.vertices, scene.triangles; color=scene.mesh_color)
end

Makie.plottype(::GeometryLike) = Plot_Geometry
end
