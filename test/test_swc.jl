@testset "test_swc.jl" begin
    swc_text = """
    # source: example
    # coordinate: micrometers

    1 1 0.0 1.0 2.0 3.0 -1
    2 3 1e0 2e0 3e0 0.5 1
    3 4 2.0 3.0 4.0 0.25 1
    # end of example
    """

    swc = mr.read_swc(IOBuffer(swc_text))
    @test swc isa mr.SWCFile
    @test swc.header == ["# source: example", "# coordinate: micrometers"]
    @test swc.footer == ["# end of example"]
    @test length(swc.nodes) == 3
    @test swc.nodes[1] == mr.SWCNode(1, 1, [0., 1., 2.], 3., -1)
    @test swc.nodes[2].position == [1., 2., 3.]
    @test swc.nodes[2].radius == 0.5
    @test swc.nodes[3].parent_id == 1

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
        loaded = mr.read_swc(filename)
        @test loaded.header == swc.header
        @test loaded.nodes == swc.nodes
        @test loaded.footer == swc.footer
    finally
        rm(filename; force=true)
    end
end
