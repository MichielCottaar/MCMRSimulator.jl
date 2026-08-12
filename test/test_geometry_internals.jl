const GI = mr.Geometries.InternalNew
const PhysicalGeometry = GI.PhysicalGeometries.PhysicalGeometry
const GeometryVector = GI.PhysicalGeometries.Groups.GeometryVector
const GeometryTuple = GI.PhysicalGeometries.Groups.GeometryTuple

@testset "Geometry internals" begin
    struct TestGeometry{N} <: PhysicalGeometry{N}
        value::Int
    end

    geometries = [TestGeometry{3}(1), TestGeometry{3}(2)]
    vector = GeometryVector{3}(geometries)
    @test vector isa GI.PhysicalGeometries.Groups.GroupGeometry{3}
    @test size(vector) == (2,)
    @test axes(vector) == (Base.OneTo(2),)
    @test length(vector) == 2
    @test firstindex(vector) == 1
    @test lastindex(vector) == 2
    @test vector[1].value == 1
    @test [g.value for g in vector] == [1, 2]
    @test [g.value for g in vector[1:2]] == [1, 2]

    vector[2] = TestGeometry{3}(3)
    @test vector[2].value == 3

    tuple = GeometryTuple{3}((TestGeometry{3}(4), TestGeometry{3}(5)))
    @test tuple isa GI.PhysicalGeometries.Groups.GroupGeometry{3}
    @test length(tuple) == 2
    @test firstindex(tuple) == 1
    @test lastindex(tuple) == 2
    @test tuple[1].value == 4
    @test [g.value for g in tuple] == [4, 5]
    @test Tuple(tuple) == (TestGeometry{3}(4), TestGeometry{3}(5))

    @test length(GeometryVector{3}(TestGeometry{3}[])) == 0
    @test length(GeometryTuple{3}(())) == 0
    @test_throws MethodError GeometryVector{3}([TestGeometry{2}(1)])
end
