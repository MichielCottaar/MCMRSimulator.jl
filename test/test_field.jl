@testset "Generate and apply microstructural fields" begin
    @testset "Simulate empty environment/sequence" begin
        spin = evolve_to_time(Spin(), Microstructure(), 0., 1.)
        @test position(spin) == SA_F64[0., 0., 0.]
        @test phase(spin) == 0.
    end
    @testset "Test constant off-resonance field" begin
        spin = evolve_to_time(Spin(), Microstructure(off_resonance=field(2.)), 0., 0.3, 3.)
        @test position(spin) == SA_F64[0., 0., 0.]
        @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
    end
    @testset "Test gradient off-resonance field" begin
        micro = Microstructure(off_resonance=field(SA_F64[1.5, 0., 0.], 2.))
        spin = evolve_to_time(Spin(), micro, 0., 0.3)
        @test position(spin) == SA_F64[0., 0., 0.]
        @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
        # Move spin and evolve a bit further in time
        spin = evolve_to_time(Spin(SA_F64[3., 0., 0.], spin.orientation), micro, 0.3, 0.5, 3.)
        @test position(spin) == SA_F64[3., 0., 0.]
        @test phase(spin) ≈ norm_angle(rad2deg((0.6 + (2. + 1.5 * 3.) * 0.2) * 3 * gyromagnetic_ratio))
    end
    @testset "Fields with different types" begin
        pos = SA_F64[1., 0., 0.]
        @test isa(field()(pos), Float64)
        @test isa(field(Int)(pos), Int)
        @test isa(field(0)(pos), Int)
        @test isa(field(0), MRSimulator.ZeroField{Int})
        @test isa(field(2)(pos), Int)
        @test isa(field(2), MRSimulator.ConstantField{Int})
        @test isa(field([0, 0, 0], 0), MRSimulator.ZeroField{Int})
        @test isa(field([0, 0, 0], 0)(pos), Int)
        @test isa(field([0, 0, 0], 2), MRSimulator.ConstantField{Int})
        @test isa(field([0, 0, 0], 2)(pos), Int)
        @test isa(field([1, 0, 0], 2), MRSimulator.GradientField{Int})
        @test isa(field([1, 0, 0], 2)(pos), Float64)
    end
end
