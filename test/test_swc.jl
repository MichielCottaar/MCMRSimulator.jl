@testset "test_swc.jl" begin
    @testset "reading an SWC file" begin
        swc_text = """
        # source: example
        # coordinate: micrometers

        1 1 0.0 1.0 2.0 3.0 -1
        2 3 1e0 2e0 3e0 0.5 1
        3 4 2.0 3.0 4.0 0.25 1
        # end of example
        """

        swc = mr.read_swc_raw(IOBuffer(swc_text))
        @test swc isa mr.SWCFile
        @test swc.header == ["# source: example", "# coordinate: micrometers"]
        @test swc.footer == ["# end of example"]
        @test length(swc.nodes) == 3
        @test swc.nodes[1] == mr.SWCNode(1, 1, [0., 1., 2.], 3., -1)
        @test swc.nodes[2].position == [1., 2., 3.]
        @test swc.nodes[2].radius == 0.5
        @test swc.nodes[3].parent_id == 1

        finite_cylinders = mr.read_geometry(IOBuffer(swc_text); format=:swc)
        @test finite_cylinders isa mr.FiniteCylinders
        @test length(finite_cylinders) == 3
        @test finite_cylinders[1].position == [0., 1., 2.]
        @test finite_cylinders[2].radius == 0.5
        @test finite_cylinders[1].connected_to == 0
        @test finite_cylinders[2].connected_to == 1
        @test mr.read_geometry(IOBuffer(swc_text)) isa mr.FiniteCylinders
        spheres_only = mr.read_geometry(IOBuffer(swc_text); format=:swc, swc_as_spheres=true)
        @test spheres_only isa mr.Spheres
        @test length(spheres_only) == 3

        ply_text = """
        ply
        format ascii 1.0
        element vertex 3
        property float x
        property float y
        property float z
        element face 1
        property list uchar int vertex_indices
        end_header
        0 0 0
        1 0 0
        0 1 0
        3 0 1 2
        """
        mesh = mr.read_geometry(IOBuffer(ply_text); format=:ply)
        @test mesh isa mr.Mesh
        @test length(mesh.vertices.value) == 3
        @test length(mesh.triangles) == 1
        @test mr.read_geometry(IOBuffer(ply_text)) isa mr.Mesh

        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 1 0\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 1 -1\n2 3 0 0 0\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 1 -1\n2 3 0 0 0 1 3\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 -1 -1\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 1 -1\n1 3 0 0 0 1 1\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("1 1 0 0 0 1 -1\n2 3 0 0 0 1 -1\n"))
        @test_throws ArgumentError mr.read_swc(IOBuffer("# header only\n"))

        filename = tempname()
        try
            open(filename, "w") do io
                write(io, swc_text)
            end
            loaded = mr.read_swc_raw(filename)
            @test loaded.header == swc.header
            @test loaded.nodes == swc.nodes
            @test loaded.footer == swc.footer
        finally
            rm(filename; force=true)
        end
    end

    @testset "reading liminal geometry files" begin
        mktempdir() do directory
            sphere_file = joinpath(directory, "sphere.json")
            swc_file = joinpath(directory, "cell.swc")
            liminal_file = joinpath(directory, "cells.txt")
            mr.write_geometry(sphere_file, mr.Spheres(radius=1.0))
            open(swc_file, "w") do io
                write(io, "1 1 0.0 0.0 0.0 1.0 -1\n")
            end
            open(liminal_file, "w") do io
                write(io, "liminal 0.2\n0.25 sphere.json\n0.75 cell.swc\n")
            end

            geometry = mr.read_geometry(liminal_file)
            @test geometry isa mr.LiminalGeometry
            @test geometry.extracellular_fraction == 0.2
            @test first.(geometry.geometries) == [0.25, 0.75]
            @test geometry.geometries[1][2] isa mr.Spheres
            @test geometry.geometries[2][2] isa mr.FiniteCylinders

            malformed_files = [
                "0.2\n1 sphere.json\n",
                "liminal nope\n1 sphere.json\n",
                "liminal 0.2\n1\n",
                "liminal 0.2\n1 missing.json\n",
            ]
            for (index, contents) in enumerate(malformed_files)
                malformed = joinpath(directory, "malformed_$index.txt")
                open(malformed, "w") do io
                    write(io, contents)
                end
                @test_throws ArgumentError mr.read_geometry(malformed)
            end
        end
    end

    @testset "reading connected SWC geometry" begin
        connected = mr.read_swc(joinpath(@__DIR__, "geometries", "cylinder.swc"), R2_inside=0.1)

        @test connected isa mr.FiniteCylinders
        @test count(==(0), connected.connected_to.value) == 1
        @test count(!=(0), connected.connected_to.value) == length(connected) - 1

        spheres = mr.read_swc(
            joinpath(@__DIR__, "geometries", "cylinder.swc"),
            R2_inside=0.1,
            swc_as_spheres=true,
        )

        seq = mr.read_pulseq(joinpath(@__DIR__, "pulseq", "gradient_echo_TE_20.seq"))

        sim = mr.Simulation(seq, geometry=spheres, diffusivity=0.5)
        snap = mr.readout(zeros(3, 300) .+ [60, 20, 20], sim, return_snapshot=true)
        @test snap.time ≈ 20.
        @test all(mr.isinside(spheres, snap) .> 0)
        @test all(mr.transverse.(snap) .≈ exp(-0.1 * 20.))

        
        snap_out = mr.readout(zeros(3, 300) .+ [58, 20, 20], sim, return_snapshot=true)
        @test snap_out.time ≈ 20.
        @test all(mr.isinside(spheres, snap_out) .== 0)
        @test all(mr.transverse.(snap_out) .≈ 1.)

    end

    @testset "diffusion through connected SWC branches" begin
        swc_text = """
        # id type x y z radius parent
        1 1 -20 0 0 10 -1
        2 3 20 0 0 10 1
        3 3 20 40 0 10 2
        4 3 20 -40 0 10 2
        """
        geometry = mr.read_swc(IOBuffer(swc_text))
        simulation = mr.Simulation([], geometry=geometry, diffusivity=3., timestep=1.)
        initial = mr.Snapshot(fill([0., 0., 0.], 1000))

        @test all(mr.isinside(geometry, initial) .> 0)

        Random.seed!(1234)
        final = mr.readout(initial, simulation, 100, return_snapshot=true)
        final_positions = mr.position.(final)
        @test any(position[2] > 15 for position in final_positions)
        @test any(position[2] < -15 for position in final_positions)
        @test all(mr.isinside(geometry, final) .> 0)
    end
end
