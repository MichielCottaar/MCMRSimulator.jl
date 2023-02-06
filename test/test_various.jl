
@test length(detect_ambiguities(mr)) == 0
@testset "Simple relaxation" begin
    orient = mr.Spin(transverse=1., longitudinal=0.).orientations[1]
    @testset "R2 relaxation" begin
        mr.relax!(orient, 0.3, mr.GlobalProperties(R2=2).mri)
        @test mr.transverse(orient) ≈ exp(-0.6)
    end
    @testset "R1 relaxation" begin
        mr.relax!(orient, 0.3, mr.GlobalProperties(R1=2).mri)
        @test mr.longitudinal(orient) ≈ 1 - exp(-0.6)
    end
end
@testset "Spin conversions" begin
    vec = SA[1, 0, 0]
    @test mr.SpinOrientation(vec).longitudinal ≈ 0.
    @test mr.SpinOrientation(vec).transverse ≈ 1.
    @test mr.phase(mr.SpinOrientation(vec)) ≈ 0.
    @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec

    vec = SA[0, 1, 1]
    @test mr.SpinOrientation(vec).longitudinal ≈ 1.
    @test mr.SpinOrientation(vec).transverse ≈ 1.
    @test mr.phase(mr.SpinOrientation(vec)) ≈ 90.
    @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec

    vec = SA[1, 1, 2]
    @test mr.SpinOrientation(vec).longitudinal ≈ 2.
    @test mr.SpinOrientation(vec).transverse ≈ sqrt(2)
    @test mr.phase(mr.SpinOrientation(vec)) ≈ 45.
    @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec
end
@testset "Apply Sequence components" begin
    @testset "0 degree flip angle pulses should do nothing" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
            pulse = mr.InstantRFPulse(0., 0., pulse_phase)
            for spin_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
                spin = mr.Spin(phase=spin_phase, transverse=1.)
                mr.apply!(pulse, spin)
                @test mr.phase(spin) ≈ spin_phase
                @test mr.longitudinal(spin) ≈ 1.
                @test mr.transverse(spin) ≈ 1.
            end
        end
    end
    @testset "180 degree pulses should flip longitudinal" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
            pulse = mr.InstantRFPulse(0., 180., pulse_phase)
            spin = mr.Spin()
            @test mr.longitudinal(spin) == 1.
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ -1.
        end
    end
    @testset "90 degree pulses should eliminate longitudinal" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
            pulse = mr.InstantRFPulse(0., 90., pulse_phase)
            spin = mr.Spin()
            @test mr.longitudinal(spin) == 1.
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
        end
    end
    @testset "Spins with same phase as pulse are unaffected by pulse" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
            for flip_angle in (10, 90, 120, 180)
                pulse = mr.InstantRFPulse(0., flip_angle, pulse_phase)

                spin = mr.Spin(longitudinal=0., transverse=1., phase=pulse_phase)
                mr.apply!(pulse, spin)
                @test mr.longitudinal(spin) ≈ zero(Float) atol=1e-7
                @test mr.transverse(spin) ≈ one(Float)
                @test mr.phase(spin) ≈ Float(pulse_phase)
            end
        end
    end
    @testset "180 pulses flips phase around axis" begin
        for spin_phase in (0., 22., 30., 80.)
            pulse = mr.InstantRFPulse(0., 180, 0.)

            spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ -spin_phase
        end
        for pulse_phase in (0., 22., 30., 80.)
            pulse = mr.InstantRFPulse(0., 180, pulse_phase)

            spin = mr.Spin(longitudinal=0., transverse=1., phase=0.)
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ 2 * pulse_phase
        end
    end
    @testset "90 pulses flips longitudinal spin into transverse plane" begin
        for pulse_phase in (0., 22., 30., 80.)
            pulse = mr.InstantRFPulse(0., 90, pulse_phase)
            spin = mr.Spin()
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ pulse_phase + 90
        end
    end
    @testset "90 pulses flips transverse spin into longitudinal plane" begin
        for pulse_phase in (0., 22., 30., 80.)
            pulse = mr.InstantRFPulse(0., 90, pulse_phase)
            spin_phase = (pulse_phase + 90)
            spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
            mr.apply!(pulse, spin)
            @test mr.longitudinal(spin) ≈ -1.
            @test mr.transverse(spin) ≈ 0. atol=1e-7
        end
    end
    @testset "Gradient should do nothing at origin" begin
        spin = mr.Spin(position=SA[0, 2, 2], transverse=1., phase=90.)
        @test mr.phase(spin) ≈ Float(90.)
        mr.apply!(mr.InstantGradient(qvec=SA[4, 0, 0]), spin)
        @test mr.phase(spin) ≈ Float(90.)
    end
    @testset "Test instant gradient effect away from origin" begin
        spin = mr.Spin(position=SA[2, 2, 2], transverse=1., phase=90.)
        @test mr.phase(spin) ≈ Float(90.)
        mr.apply!(mr.InstantGradient(qvec=SA[0.01, 0, 0]), spin)
        @test mr.phase(spin) ≈ Float(90. + 0.02 * 360)
    end
end
@testset "Random generator number control" begin
    @testset "FixedXoshiro interface" begin
        rng = mr.FixedXoshiro()
        copy!(Random.TaskLocalRNG(), rng)
        a = randn(2)
        copy!(Random.TaskLocalRNG(), rng)
        @test randn(2) == a
    end
    @testset "FixedXoshiro predictability" begin
        Random.seed!(1234)
        rng = mr.FixedXoshiro()
        Random.seed!(1234)
        rng2 = mr.FixedXoshiro()
        @test rng == rng2
    end
    @testset "Reproducible evolution" begin
        spin = mr.Spin()
        env = mr.Simulation([], diffusivity=3.)
        t1 = [s[1] for s in mr.trajectory(spin, env, 1:5)]
        t2 = [s[1] for s in mr.trajectory(spin, env, 1:5)]
        @test all(@. mr.position(t1) == mr.position(t2))
        get_rng(spin) = spin.rng
        @test all(@. get_rng(t1) == get_rng(t2))
    end
end
@testset "Bounding boxes" begin
    # single obstruction
    @test mr.BoundingBox(mr.Cylinder(1)) == mr.BoundingBox([-1, -1], [1, 1])

    # shifted cylinder
    @test mr.BoundingBox(mr.cylinders(1., positions=[2., 2.])) == mr.BoundingBox([1., 1.], [3, 3])

    # repeated obstructions
    @test mr.BoundingBox(mr.cylinders(1., repeats=[1., 3.])) == mr.BoundingBox([-1, -1], [1, 1])

    # shifted spheres
    @test mr.BoundingBox(mr.spheres(1., positions=[[1, 0, 0], [0, 1, 0]])) == mr.BoundingBox([-1, -1, -1.], [2., 2., 1.])
end

@testset "Test readout formats" begin
    seq = mr.Sequence(TR=100, pulses=[mr.Readout(10.), mr.Readout(30.)])
    sim_empty = mr.Simulation([])
    sim_flat = mr.Simulation(seq)
    sim_single = mr.Simulation([seq])
    sim_double = mr.Simulation([seq, seq])

    @testset "Test signal function output" begin
        times = 1:0.1:10
        Nt = length(times)
        for (sim, shape) in [
            (sim_empty, (Nt, 0)),
            (sim_flat, (Nt, )),
            (sim_single, (Nt, 1)),
            (sim_double, (Nt, 2)),
        ]
            @test size(mr.signal(100, sim, times)) == shape
            @test all(mr.longitudinal.(mr.signal(100, sim, times)) .≈ 100.)
        end
    end
    @testset "Test readout function output" begin
        @test length(mr.readout(100, sim_empty)) == 0
        @test length(mr.readout(100, sim_flat)) == 2
        for (time, snapshot) in zip((10., 30.), mr.readout(100, sim_flat))
            @test length(snapshot) == 100
            @test time == mr.get_time(snapshot)
            @test mr.longitudinal(snapshot) ≈ 100.
        end
        @test length(mr.readout(100, sim_single)) == 1
        @test length(mr.readout(100, sim_double)) == 2
        for sim in (sim_single, sim_double)
            @test length(sim.sequences) == length(mr.readout(100, sim))
            for r in mr.readout(100, sim)
                @test length(r) == 2
                for (time, snapshot) in zip((10., 30.), r)
                    @test length(snapshot) == 100
                    @test time == mr.get_time(snapshot)
                    @test mr.longitudinal(snapshot) ≈ 100.
                end
            end
        end
    end
end

@testset "Test simulation pretty printing" begin
    sim = mr.Simulation([mr.dwi(bval=1), mr.dwi(bval=2)], geometry=mr.spheres([1, 2.], repeats=[5, 5, 5]), R1=0.1, MT_fraction=0.3)
    @test repr(sim, context=:compact => true) == "Simulation(2 sequences, Geometry(2 repeating spheres, ), D=0.0um^2/ms, GlobalProperties(T1=10.0ms, MT_fraction=0.3, ))" 
    @test repr(sim, context=:compact => false) == "Simulation(Geometry(2 repeating spheres, ), D=0.0um^2/ms, GlobalProperties(T1=10.0ms, MT_fraction=0.3, )):
2 sequences:
Sequence (TR=2000.0ms):
    - InstantRFPulse: t=0.0ms, θ=90.0°, ϕ=-90.0°;
    - InstantRFPulse: t=40.0ms, θ=180.0°, ϕ=0.0°;
    - Readout at 80.0ms
Sequence (TR=2000.0ms):
    - InstantRFPulse: t=0.0ms, θ=90.0°, ϕ=-90.0°;
    - InstantRFPulse: t=40.0ms, θ=180.0°, ϕ=0.0°;
    - Readout at 80.0ms
" 
end