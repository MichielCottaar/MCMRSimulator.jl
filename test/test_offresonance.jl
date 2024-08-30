
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
            cylinders = mr.Cylinders(radius=1.)
            @test field(cylinders, zero(SVector{3, Float64})) == 0.
            @test field(cylinders, SVector{3, Float64}([1, 1, 1])) == 0.
            @test iszero(grad(cylinders))
        end
        @testset "Myelinated cylinders aligned with field produce no off-resonance" begin
            cylinders = mr.Cylinders(radius=1., position=[[0, 0], [2, 0]], g_ratio=0.8)
            @test field(cylinders, zero(SVector{3, Float64})) == 0.
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) == 0.
            @test iszero(grad(cylinders))
        end
        @testset "Myelinated cylinder perpendicular to field with only isotropic susceptibility" begin
            cylinders = mr.Cylinders(radius=1., rotation=:x, susceptibility_aniso=0., susceptibility_iso=1., g_ratio=0.8, grid_resolution=Inf)
            @test field(cylinders, zero(SVector{3, Float64})) ≈ 0.
            @test field(cylinders, [0., 0.1, 0.]) ≈ 0.
            outer_field = 1//2 * (1 - 0.8^2) / (1 + 0.8)^2
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(cylinders) ≈ outer_field * 4
        end
        @testset "Myelinated cylinder perpendicular to field with only anisotropic susceptibility" begin
            cylinders = mr.Cylinders(radius=1., rotation=:x, susceptibility_aniso=1., susceptibility_iso=0., g_ratio=0.8, grid_resolution=Inf)
            @test field(cylinders, zero(SVector{3, Float64})) ≈ -0.75 * log(0.8)
            outer_field = 1//8 * (1 - 0.8^2) / (1 + 0.8)^2
            @test field(cylinders, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(cylinders, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(cylinders) ≈ outer_field * 4
        end
        @testset "Repeating myelinated cylinder perpendicular to field with only isotropic susceptibility" begin
            cylinders = mr.Cylinders(radius=1., rotation=:x, susceptibility_aniso=0., susceptibility_iso=1., g_ratio=0.8, repeats=(4, 4), grid_resolution=Inf)
            @test field(cylinders, zero(SVector{3, Float64})) ≈ 0. atol=1e-12
            @test field(cylinders, [0., 0., 0.2]) > 1e-6
            @test field(cylinders, [0., 0.2, 0.]) < -1e-6
            outer_field = 1//2 * (1 - 0.8^2) / (1 + 0.8)^2
            @test field(cylinders, [0, 0, 2]) ≈ outer_field * 2 rtol=0.2
            @test field(cylinders, [0, 2, 0]) ≈ -outer_field * 2 rtol=0.2
            @test grad(cylinders) ≈ outer_field * 4
        end
    end
    @testset "Annulus off-resonance field" begin
        @testset "Unmyelinated annulus produces no field" begin
            annuli = mr.Annuli(inner=0.7, outer=1.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([1, 1, 1])) == 0.
            @test iszero(grad(annuli))
        end
        @testset "Myelinated annuli aligned with field only produce off-resonance within the myelin" begin
            annuli = mr.Annuli(inner=0.7, outer=1., position=[[0, 0], [2, 0]], myelin=true, susceptibility_iso=1., susceptibility_aniso=0.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0.8, 0])) ≈ 1//3
            @test field(annuli, SVector{3, Float64}([1.2, 0, 0])) ≈ 1//3

            annuli = mr.Annuli(inner=0.7, outer=1., position=[[0, 0], [2, 0]], myelin=true, susceptibility_iso=0., susceptibility_aniso=1.)
            @test field(annuli, zero(SVector{3, Float64})) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) == 0.
            @test field(annuli, SVector{3, Float64}([0, 0.8, 0])) ≈ -1//6
            @test field(annuli, SVector{3, Float64}([1.2, 0, 0])) ≈ -1//6
        end
        @testset "Myelinated annulus perpendicular to field with only isotropic susceptibility" begin
            annuli = mr.Annuli(inner=0.7, outer=1., rotation=:x, susceptibility_aniso=0., susceptibility_iso=1., myelin=true, grid_resolution=Inf)
            @test field(annuli, zero(SVector{3, Float64})) ≈ 0.
            outer_field = 1//2 * (1 - 0.7^2) / 4
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(annuli) ≈ outer_field * 4
        end
        @testset "Myelinated annulus perpendicular to field with only anisotropic susceptibility" begin
            annuli = mr.Annuli(inner=0.7, outer=1., rotation=:x, susceptibility_aniso=1., susceptibility_iso=0., myelin=true, grid_resolution=Inf)
            @test field(annuli, zero(SVector{3, Float64})) ≈ -0.75 * log(0.7)
            outer_field = 1//8 * (1 - 0.7^2) / 4
            @test field(annuli, SVector{3, Float64}([0, 0, 2])) ≈ outer_field
            @test field(annuli, SVector{3, Float64}([0, 2, 0])) ≈ -outer_field
            @test grad(annuli) ≈ outer_field * 4
        end
    end
    @testset "mesh off-resonance field" begin
        for aniso in (0, 1e-12) # test isotropic and anisotropic version of the code
            @testset "Single right triangle test" begin
                for t in ([1, 2, 3], [2, 1, 3])
                    mesh = mr.Mesh(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0]], triangles=[t], myelin=true, susceptibility_iso=1., susceptibility_aniso=aniso)
                    f(z) = atan(sqrt(2 + z^2)/ z)
                    ref_field = (f(0.99) - f(1.01)) / (0.02 * 4π)
                    @test field(mesh, [0, 0, 1]) ≈ ref_field rtol=0.1
                    @test field(mesh, [0, 0, -1]) ≈ ref_field rtol=0.1
                end
            end
            @testset "Points very close to mesh" begin
                for t in ([1, 2, 3], [2, 1, 3])
                    mesh = mr.Mesh(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0]], triangles=[t], myelin=true, susceptibility_iso=1., susceptibility_aniso=aniso)
                    for y in (0.2, 0.3)
                        ref = field(mesh, [0.8, y, 1e-4])
                        @test field(mesh, [0.8, y, 1e-5]) ≈ ref rtol=0.01
                        @test field(mesh, [0.8, y, 1e-3]) ≈ ref rtol=0.01
                        @test field(mesh, [0.8, y, -1e-5]) ≈ ref rtol=0.01
                        @test field(mesh, [0.8, y, -1e-3]) ≈ ref rtol=0.01
                    end
                end
            end
            @testset "Points very far above mesh" begin
                for t in ([1, 2, 3], [2, 1, 3])
                    for v in [[1, 1, 0], [0, 1, 0]]
                        mesh = mr.Mesh(vertices=[[0, 0, 0], [1, 0, 0], v], triangles=[t], myelin=true, susceptibility_iso=1., susceptibility_aniso=aniso, lorentz_radius=1e7, grid_resolution=Inf)
                        @test field(mesh, [1e-5, 1e-5, 20.]) ≈ 1 / (20^3 * 4π) rtol=0.1
                        @test field(mesh, [1e-5, 1e-5, -20.]) ≈ 1 / (20^3 * 4π) rtol=0.1
                    end
                end
            end
            @testset "Points very to side of mesh" begin
                for t in ([1, 2, 3], [2, 1, 3])
                    for v in [[1, 0, 1], [0, 0, 1]]
                        mesh = mr.Mesh(vertices=[[0, 0, 0], [1, 0, 0], v], triangles=[t], myelin=true, susceptibility_iso=1., susceptibility_aniso=aniso, lorentz_radius=1e7, grid_resolution=Inf)
                        @test field(mesh, [1e-3, 1e-3, 20.]) ≈ 1 / 4π / 20^3 rtol=0.1
                        @test field(mesh, [1e-3, 1e-3, -20.]) ≈ 1 / 4π / 20^3 rtol=0.1
                    end
                end
            end
        end
    end
end