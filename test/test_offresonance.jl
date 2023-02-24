
@testset "test_offresonance.jl" begin
    @testset "Cylinder off-resonance field" begin
        @testset "Unmyelinated cylinder produces no field" begin
            cylinders = mr.cylinders(1.)
            @test mr.off_resonance(cylinders, zero(mr.PosVector)) == 0.
            @test mr.off_resonance(cylinders, mr.PosVector([1, 1, 1])) == 0.
        end
        @testset "Myelinated cylinders aligned with field produce no off-resonance" begin
            cylinders = mr.cylinders(1., positions=[[0, 0], [2, 0]], g_ratio=0.8)
            @test mr.off_resonance(cylinders, zero(mr.PosVector)) == 0.
            @test mr.off_resonance(cylinders, mr.PosVector([0, 0, 2])) == 0.
            @test mr.off_resonance(cylinders, mr.PosVector([0, 2, 0])) == 0.
        end
        @testset "Myelinated cylinder perpendicular to field with only isotropic susceptibility" begin
            cylinders = mr.cylinders(1., rotation=:x, chi_A=0., chi_I=1., g_ratio=0.8)
            @test mr.off_resonance(cylinders, zero(mr.PosVector)) ≈ 0.
            outer_field = 1//2 * (1 - 0.8^2) / (1 + 0.8)^2
            @test mr.off_resonance(cylinders, mr.PosVector([0, 0, 2])) ≈ outer_field
            @test mr.off_resonance(cylinders, mr.PosVector([0, 2, 0])) ≈ -outer_field
        end
        @testset "Myelinated cylinder perpendicular to field with only anisotropic susceptibility" begin
            cylinders = mr.cylinders(1., rotation=:x, chi_A=1., chi_I=0., g_ratio=0.8)
            @test mr.off_resonance(cylinders, zero(mr.PosVector)) ≈ -0.75 * log(0.8)
            outer_field = 1//8 * (1 - 0.8^2) / (1 + 0.8)^2
            @test mr.off_resonance(cylinders, mr.PosVector([0, 0, 2])) ≈ outer_field
            @test mr.off_resonance(cylinders, mr.PosVector([0, 2, 0])) ≈ -outer_field
        end
    end
end