@testset "Generate and apply microstructural fields" begin
    function mr.evolve_to_time(
        spin::mr.Spin{N}, current_time::Float, new_time::Float,
        micro::mr.Microstructure, timestep::Float, B0::Float=Float(3.)
    ) where {N}
        empty_pulses = SVector{N}([(mr.InstantRFPulse(0, 0, 0), mr.InstantRFPulse(0, 0, 0))])
        spin2 = copy(spin)
        mr.evolve_to_time!(spin2, mr.Simulation([mr.Sequence(TR=new_time, B0=B0) for _ in 1:N], micro, mr.TimeController(timestep)), empty_pulses, current_time, new_time)
        spin2
    end

    for step_size in (1., 0.3, 0.123)
        step_size = Float(step_size)
        @testset "Simulate empty environment/sequence with dt=$step_size" begin
            orient = mr.evolve_to_time(mr.Spin(transverse=1.), zero(Float), one(Float), mr.Microstructure(), step_size).orientations[1]
            @test mr.phase(orient) == 0.
        end
        @testset "Test constant off-resonance field with dt=$step_size" begin
            orient = mr.evolve_to_time(mr.Spin(transverse=1.), zero(Float), Float(0.3), mr.Microstructure(off_resonance=mr.field(2.)), step_size).orientations[1]
            @test mr.phase(orient) ≈ mr.norm_angle(0.3 * 2. * 360)
        end
        @testset "Test gradient off-resonance field with dt=$step_size" begin
            micro = mr.Microstructure(off_resonance=mr.field(SA[1.5, 0., 0.], 2.))
            spin = mr.evolve_to_time(mr.Spin(transverse=1.), zero(Float), Float(0.3), micro, step_size)
            orient = spin.orientations[1]
            @test mr.phase(orient) ≈ mr.norm_angle(0.6 * 360)
            # Move spin and evolve a bit further in time
            orient = mr.evolve_to_time(mr.Spin(transverse=1., position=SA[3., 0., 0.]), zero(Float), Float(0.2), micro, step_size).orientations[1]
            @test mr.phase(orient) ≈ mr.norm_angle((2. + 1.5 * 3.) * 0.2 * 360)
        end
    end
    @testset "Fields with different types" begin
        pos = SA[1., 0., 0.]
        @test isa(mr.field()(pos), Float)
        @test isa(mr.field(0)(pos), Float)
        @test isa(mr.field(0), mr.ZeroField)
        @test isa(mr.field(2)(pos), Float)
        @test isa(mr.field(2), mr.ConstantField)
        @test isa(mr.field([0, 0, 0], 0), mr.ZeroField)
        @test isa(mr.field([0, 0, 0], 0)(pos), Float)
        @test isa(mr.field([0, 0, 0], 2), mr.ConstantField)
        @test isa(mr.field([0, 0, 0], 2)(pos), Float)
        @test isa(mr.field([1, 0, 0], 2), mr.GradientField)
        @test isa(mr.field([1, 0, 0], 2)(pos), Float)
    end
end
