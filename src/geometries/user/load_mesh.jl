module LoadMesh
import PlyIO
import ..Obstructions: mesh
"""
    ply_from_mesh(file)

Loads a [`Mesh`](@ref) from a PLY file.
PLY stands for Polygon File Format (http://paulbourke.net/dataformats/ply/).
PLY IO is handled by PlyIO.jl (https://github.com/JuliaGeometry/PlyIO.jl).
"""
function ply_from_mesh(ply_file)
    ply = PlyIO.load_ply(ply_file)
    vertices = zip(
        ply["vertex"]["x"],
        ply["vertex"]["y"],
        ply["vertex"]["z"],
    )
    indices = [v .+ 1 for v in ply["face"]["vertex_indices"]]
    return mesh(vertices=vertices, triangles=indices)
end


"""
    load_mesh(file)

Loads a [`Mesh`](@ref) from a file.

Currently only PLY files are supported (see [`ply_from_mesh`](@ref))
"""
function load_mesh(io::IO)
    mark(io)
    firstline = readline(io)
    reset(io)
    if firstline == "ply"
        return ply_from_mesh(io)
    else
        throw(ErrorException("Could not recognise input mesh format!"))
    end
end


function load_mesh(file_name::AbstractString)
    open(file_name, "r") do fid
        load_mesh(fid)
    end
end

end