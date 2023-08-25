
@testset "test_offresonance.jl" begin
    function field(geometry, position)
        susc = mr.fix_susceptibility(geometry)
        return mr.Geometries.Internal.susceptibility_off_resonance(susc, SVector{3, Float64}(position))
    end
    function grad(geometry)
        susc = mr.fix_susceptibility(geometry)
        return mr.Geometries.Internal.off_resonance_gradient(susc, 1/mr.gyromagnetic_ratio) * 1e6
    end
    @testset "Cylinder off-resonance field" begin
        @testset "Unmyelinated cylinder produces no field" begin
            cylinders = mr.cylinders(radius=1.)
            @test field(cylinders, zero(SVector{3, Float64})) == 0.
            @test field(cylinders, SVector{3, Float64}([1, 1, 1])) == 0.
            @test iszero(grad(cylinders))
        end
        @testset "Myelinated cylinders aligned with field produce no off-resonance" begin
            cylinders = mr.cylinders(radius=1., position=[[0, 0], [2, 0]], g_ratio=0.8)
            @test field(cylinders, zero(SVector{3, Float64})) == 0.
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) == 0.
            @test iszero(grad(cylinders))
        end
        @testset "Myelinated cylinder perpendicular to field with only isotropic susceptibility" begin
            cylinders = mr.cylinders(radius=1., rotation=:x, susceptibility_aniso=0., susceptibility_iso=1., g_ratio=0.8)
            @test field(cylinders, zero(SVector{3, Float64})) ≈ 0.
            outer_field = 1//2 * (1 - 0.8^2) / (1 + 0.8)^2
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(cylinders) ≈ outer_field * 4
        end
        @testset "Myelinated cylinder perpendicular to field with only anisotropic susceptibility" begin
            cylinders = mr.cylinders(radius=1., rotation=:x, susceptibility_aniso=1., susceptibility_iso=0., g_ratio=0.8)
            @test field(cylinders, zero(SVector{3, Float64})) ≈ -0.75 * log(0.8)
            outer_field = 1//8 * (1 - 0.8^2) / (1 + 0.8)^2
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(cylinders) ≈ outer_field * 4
        end
    end
    @testset "Annulus off-resonance field" begin
        @testset "Unmyelinated annulus produces no field" begin
            annuli = mr.annuli(inner=0.7, outer=1.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([1, 1, 1])) == 0.
            @test iszero(grad(annuli))
        end
        @testset "Myelinated annuli aligned with field only produce off-resonance within the myelin" begin
            annuli = mr.annuli(inner=0.7, outer=1., position=[[0, 0], [2, 0]], myelin=true, susceptibility_iso=1., susceptibility_aniso=0.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0.8, 0])) ≈ 1//3
            @test field(annuli, SVector{3, Float64}([1.2, 0, 0])) ≈ 1//3

            annuli = mr.annuli(inner=0.7, outer=1., position=[[0, 0], [2, 0]], myelin=true, susceptibility_iso=0., susceptibility_aniso=1.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0.8, 0])) ≈ -1//6
            @test field(annuli, SVector{3, Float64}([1.2, 0, 0])) ≈ -1//6
        end
        @testset "Myelinated annulus perpendicular to field with only isotropic susceptibility" begin
            annuli = mr.annuli(inner=0.7, outer=1., rotation=:x, susceptibility_aniso=0., susceptibility_iso=1., myelin=true)
            @test field(annuli, zero(SVector{3, Float64})) ≈ 0.
            outer_field = 1//2 * (1 - 0.7^2) / 4
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(annuli) ≈ outer_field * 4
        end
        @testset "Myelinated annulus perpendicular to field with only anisotropic susceptibility" begin
            annuli = mr.annuli(inner=0.7, outer=1., rotation=:x, susceptibility_aniso=1., susceptibility_iso=0., myelin=true)
            @test field(annuli, zero(SVector{3, Float64})) ≈ -0.75 * log(0.7)
            outer_field = 1//8 * (1 - 0.7^2) / 4
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(annuli) ≈ outer_field * 4
        end
    end
    @testset "mesh off-resonance field" begin
        @testset "Single right triangle test" begin
            for t in ([1, 2, 3], [2, 1, 3])
                mesh = mr.mesh(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0]], triangles=[t], myelin=true, susceptibility_iso=1., susceptibility_aniso=0.)
                @test field(mesh, [0, 0, 1]) ≈ atan(sqrt(3)) - atan(1)
                @test field(mesh, [0, 0, -1]) ≈ atan(-sqrt(3)) + atan(1)
            end
        end
    end
end