@testset "test_cli.jl" begin
import MCMRSimulator.CLI: run_main_test
import JSON
import CSV 
import DataFrames: DataFrame

function in_tmpdir(f)
    curdir = pwd()
    mktempdir() do fn
        cd(fn)
        try
            f()
        finally
            cd(curdir)
        end
    end
end

@testset "Test the creation commands do not crash" begin
    @testset "mcmr geometry create" begin
        @testset "create cylinders" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create cylinders 2 test.json --radius 1,2 --position 0,1;2,3 --R1_surface 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test result["number"] == 2
                @test result["radius"] == [1., 2.]
                @test result["position"] == [[0., 1.], [2., 3.]]
                @test result["R1_surface"] == 1.5
                @test result["R1_inside"] == 0.
            end
        end
        @testset "create spheres" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create spheres 2 test.json --radius 1,2 --position 0,1,2;3,4,5 --R1_surface 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test result["number"] == 2
                @test result["radius"] == [1., 2.]
                @test result["position"] == [[0., 1., 2.], [3., 4., 5.]]
                @test result["R1_surface"] == 1.5
                @test result["R1_inside"] == 0.
            end
        end
        @testset "create walls" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create walls 2 test.json --position 0,1 --R1 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Walls"
                @test result["number"] == 2
                @test result["position"] == [0., 1.]
                @test result["R1"] == 1.5
                @test result["R2"] == 0.
            end
        end
        @testset "create annuli" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create annuli 2 test.json --inner 1,2 --outer 3 --myelin")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test result["number"] == 2
                @test result["position"] == [0., 0.]
                @test result["inner"] == [1, 2]
                @test result["outer"] == 3.
                @test result["myelin"] == true
            end
        end
        @testset "create bendy cylinder" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create bendy-cylinder 2 test.json --control-point 0,0,0;0,0.5,1 --radius 0.4,0.6 --closed 0,0,1 --repeats 1.3,2,2")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "BendyCylinder"
                @test result["number"] == 2
                @test result["position"] == [[0., 0., 0.], [0., 0.5, 1.]]
                @test result["radius"] == [0.4, 0.6]
                @test result["closed"] == [0, 0, 1]
                @test result["repeats"] == [1.3, 2., 2.]
            end
        end
        @testset "create based on file" begin
            in_tmpdir() do 
                open("positions.txt", "w") do f
                    write(f, "0 0
                    1 1
                    2 0.5")
                end
                open("radii.txt", "w") do f
                    write(f, "1.2 0.8 0.3")
                end
                open("R1.txt", "w") do f
                    write(f, "0.1")
                end
                _, err = run_main_test("geometry create cylinders 3 test.json --position positions.txt --radius radii.txt --R1_inside R1.txt")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test result["number"] == 3
                @test result["position"] == [[0., 0.], [1., 1.], [2, 0.5]]
                @test result["radius"] == [1.2, 0.8, 0.3]
                @test result["R1_inside"] == 0.1
                @test result["R1_surface"] == 0.
            end
        end
    end
    @testset "mcmr geometry create-random" begin
        @testset "random cylinders" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random cylinders 0.5 test.json --repeats 20,20")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random spheres" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random spheres 0.5 test.json --repeats 20,20,20")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random annuli" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random annuli 0.5 test.json --repeats 20,20 --g-ratio=0.8")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test all(result["outer"] .≈ 1.)
                @test all(result["inner"] .≈ 0.8)
                @test result["number"] > 3
            end
        end
    end
    @testset "mcmr sequence" begin
        @testset "gradient_echo" begin
            in_tmpdir() do
                _, err = run_main_test("sequence gradient-echo seq.json --TE 30 --TR 100")
                @test length(err) == 0
                sequence = mr.read_sequence("seq.json")
                @test length(sequence.pulses) == 0
                @test length(sequence.instants) == 2
                @test sequence.instants[1].time == 0.
                @test sequence.instants[2].time > 30.
                @test sequence.instants[1] isa mr.InstantRFPulse
                @test sequence.instants[2] isa mr.InstantGradient
                @test length(sequence.gradients) == 0
                @test sequence.TR == 100
                @test sequence.readout_times == [30.]
            end
        end
        @testset "spin_echo" begin
            in_tmpdir() do
                _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100")
                @test length(err) == 0
                sequence = mr.read_sequence("seq.json")
                @test length(sequence.pulses) == 0
                @test length(sequence.instants) == 3
                @test sequence.instants[1].time == 0.
                @test sequence.instants[2].time == 15.
                @test sequence.instants[3].time > 30.
                @test sequence.instants[1] isa mr.InstantRFPulse
                @test sequence.instants[2] isa mr.InstantRFPulse
                @test sequence.instants[3] isa mr.InstantGradient
                @test length(sequence.gradients) == 0
                @test sequence.TR == 100
                @test sequence.readout_times == [30.]
            end
        end

        @testset "Test controls of crusher gradients" begin
            @testset "Spin echo sequence crushers" begin
                in_tmpdir() do
                    _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100 --refocus-crusher-duration 2 --crusher-duration=1 --scanner=Siemens_Prisma --readout-time=1")
                    @test length(err) == 0
                    sequence = mr.read_sequence("seq.json")
                    @test length(sequence.pulses) == 0
                    @test length(sequence.instants) == 2
                    @test sequence.instants[1].time == 0.
                    @test sequence.instants[2].time == 15.
                    @test sequence.instants[1] isa mr.InstantRFPulse
                    @test sequence.instants[2] isa mr.InstantRFPulse
                    @test length(sequence.gradients) == 3
                    @test sequence.gradients[1].shape.times[1] == 13.
                    @test sequence.gradients[1].shape.times[end] == 15
                    @test sequence.gradients[2].shape.times[1] == 15.
                    @test sequence.gradients[2].shape.times[end] == 17.
                    @test sequence.gradients[3].shape.times[1] == 30.5
                    @test sequence.gradients[3].shape.times[end] == 31.5
                    for grad in sequence.gradients
                        amplitude = norm(grad.shape.amplitudes[2])
                        @test amplitude ≈ √3 * mr.max_gradient(mr.Siemens_Prisma)
                        @test amplitude / (grad.shape.times[2] - grad.shape.times[1]) ≈ √3 * mr.max_slew_rate(mr.Siemens_Prisma)
                    end
                    @test sequence.TR == 100
                    @test sequence.readout_times == [30.]
                end
            end
            @testset "Spin echo sequence crushers using q-val" begin
                @testset "Using infinite scanner" begin
                    in_tmpdir() do
                        _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100 --refocus-crusher-qval 2 --crusher-duration=4 --readout-time=1")
                        @test length(err) == 0
                        sequence = mr.read_sequence("seq.json")
                        @test length(sequence.instants) == 4
                        @test sequence.instants[1].time == 0.
                        @test 14.99 < sequence.instants[2].time < 15.
                        @test sequence.instants[3].time == 15.
                        @test 15. < sequence.instants[4].time < 15.01
                        @test sequence.instants[1] isa mr.InstantRFPulse
                        @test sequence.instants[2] isa mr.InstantGradient
                        @test sequence.instants[3] isa mr.InstantRFPulse
                        @test sequence.instants[4] isa mr.InstantGradient
                        @test length(sequence.gradients) == 1
                        grad = sequence.gradients[1]
                        @test grad.shape.times[1] == 30.5
                        @test 0 < (grad.shape.times[2] - 30.5) < 1e-3
                        @test -1e-3 < (grad.shape.times[3] - 34.5) < 0
                        @test grad.shape.times[end] == 34.5
                        @test sequence.TR == 100
                        @test sequence.readout_times == [30.]
                    end
                end
                @testset "Using finite scanner" begin
                    in_tmpdir() do
                        _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100 --refocus-crusher-qval 0.02 --crusher-duration=4 --readout-time=1 --scanner=Siemens_Prisma")
                        @test length(err) == 0
                        sequence = mr.read_sequence("seq.json")
                        @test length(sequence.instants) == 2
                        @test sequence.instants[1].time == 0.
                        @test sequence.instants[2].time == 15.
                        @test sequence.instants[1] isa mr.InstantRFPulse
                        @test sequence.instants[2] isa mr.InstantRFPulse
                        @test length(sequence.gradients) == 3
                        @test sequence.gradients[1].shape.times[end] == 15.
                        @test sequence.gradients[2].shape.times[1] == 15.
                        for grad in sequence.gradients[1:2]
                            amplitude = norm(grad.shape.amplitudes[2])
                            @test amplitude ≈ √3 * mr.max_gradient(mr.Siemens_Prisma)
                            qval = (grad.shape.times[3] - grad.shape.times[1]) * amplitude
                            @test qval ≈ 0.02
                        end
                        grad = sequence.gradients[3]
                        @test grad.shape.times[1] ≈ 30.5
                        @test grad.shape.times[end] == 34.5
                        @test sequence.TR == 100
                        @test all(sequence.readout_times .≈ [30.])
                    end
                end
            end
            @testset "gradient echo crushers" begin
                in_tmpdir() do
                    _, err = run_main_test("sequence gradient-echo seq.json --TE 30 --TR 100 --crusher-duration=1 --scanner=Siemens_Connectom --readout-time=1")
                    @test length(err) == 0
                    sequence = mr.read_sequence("seq.json")
                    @test length(sequence.pulses) == 0
                    @test length(sequence.instants) == 1
                    @test sequence.instants[1].time == 0.
                    @test sequence.instants[1] isa mr.InstantRFPulse
                    @test length(sequence.gradients) == 1
                    grad = sequence.gradients[1]
                    @test grad.shape.times[1] == 30.5
                    # maximum amplitude has not been reached
                    @test grad.shape.times[2] ≈ 31
                    @test grad.shape.times[3] ≈ 31
                    @test grad.shape.times[end] == 31.5
                    amplitude = norm(grad.shape.amplitudes[2])
                    # maximum amplitude has not been reached
                    @test amplitude < √3 * mr.max_gradient(mr.Siemens_Connectom)
                    @test amplitude / (grad.shape.times[2] - grad.shape.times[1]) ≈ √3 * mr.max_slew_rate(mr.Siemens_Connectom)
                    @test sequence.TR == 100
                    @test sequence.readout_times == [30.]
                end
            end
        end
        @testset "spin_echo with finite pulses" begin
            in_tmpdir() do
                _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100 --excitation-duration 2 --refocus-duration 4")
                @test length(err) == 0
                sequence = mr.read_sequence("seq.json")
                @test length(sequence.pulses) == 2
                @test length(sequence.instants) == 1
                @test length(sequence.gradients) == 0
                @test sequence.TR == 100
                @test sequence.readout_times == [31.]
            end
        end
        @testset "DWI" begin
            in_tmpdir() do
                _, err = run_main_test("sequence dwi seq.json --TE 30 --TR 100 --bval 2")
                @test length(err) == 0
                sequence = mr.read_sequence("seq.json")
                @test length(sequence.pulses) == 0
                @test length(sequence.instants) == 3
                @test sequence.instants[1].time == 0.
                @test sequence.instants[2].time ≈ 15.
                @test sequence.instants[3].time > 30.
                @test sequence.instants[1] isa mr.InstantRFPulse
                @test sequence.instants[2] isa mr.InstantRFPulse
                @test sequence.instants[3] isa mr.InstantGradient
                @test length(sequence.gradients) == 2
                @test sequence.TR == 100
                @test sequence.readout_times == [30.]
            end
        end
    end
end
@testset "Test full simulations" begin
    @testset "Setting R2" begin
        in_tmpdir() do
            @testset "Single global R2" begin
                _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2")
                @test length(err) == 0
                _, err = run_main_test("sequence gradient-echo ge.json --TE 30 --TR 100")
                @test length(err) == 0
                _, err = run_main_test("run spheres.json ge.json --R2 0.1 -N 100 -o global.csv")
                @test length(err) == 0
                result = DataFrame(CSV.File("global.csv"))
                @test size(result, 1) == 1
                @test result[1, :nspins] == 100
                @test result[1, :transverse] ≈ 100 * exp(-3.)
            end
            @testset "Change R2 within sphere" begin
                _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2 --R2_inside=0.2")
                @test length(err) == 0
                _, err = run_main_test("run spheres.json ge.json --R2 0.1 -N 100 -o varies.csv --subset inside --subset outside")
                @test length(err) == 0
                result = DataFrame(CSV.File("varies.csv"))
                @test size(result, 1) == 3
                @test result[1, :nspins] == 100
                @test result[2, :nspins] < 100
                @test result[3, :nspins] < 100
                @test result[2, :transverse] / result[2, :nspins] ≈ exp(-9.)
                @test result[3, :transverse] / result[3, :nspins] ≈ exp(-3.)
                @test all(result[!, :sequence] .== "ge.json")
                @test all(result[!, :sequence_index] .== 1)
            end
        end
    end
end
@testset "Setting R1" begin
    in_tmpdir() do
        _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2 --R1_inside=0.02")
        @test length(err) == 0
        _, err = run_main_test("sequence gradient-echo ge.json --TE 30 --TR 100")
        @test length(err) == 0
        _, err = run_main_test("run spheres.json ge.json --R1 0.01 -N 100 -o R1.csv --subset inside --subset outside")
        @test length(err) == 0
        result = DataFrame(CSV.File("R1.csv"))
        @test size(result, 1) == 3
        @test result[1, :nspins] == 100
        @test result[2, :nspins] < 100
        @test result[3, :nspins] < 100
        @test result[2, :longitudinal] / result[2, :nspins] ≈ 1. - exp(-0.9)
        @test result[3, :longitudinal] / result[3, :nspins] ≈ 1. - exp(-0.3)
        @test all(result[!, :sequence] .== "ge.json")
        @test all(result[!, :sequence_index] .== 1)
    end
end
@testset "Setting diffusivity" begin
    in_tmpdir() do
        _, err = run_main_test("geometry create walls 1 walls.json --repeats 1")
        @test length(err) == 0
        _, err = run_main_test("sequence dwi dwi.json --bval 0.5 --TE 80 --TR 1000")
        @test length(err) == 0
        open("bvecs", "w") do f
            write(f, "1 0
            0 1
            0 0")
        end
        _, err = run_main_test("run walls.json dwi.json -o with_diff.csv --diffusivity 3. --bvecs=bvecs")
        @test length(err) == 0
        with_diff = DataFrame(CSV.File("with_diff.csv"))
        @test size(with_diff, 1) == 2
        @test with_diff[!, :sequence] == [1, 1]
        @test with_diff[!, :bvec] == [1, 2]
        @test with_diff[1, :transverse] > exp(-1) * with_diff[1, :nspins]
        @test with_diff[2, :transverse] ≈ with_diff[2, :nspins] * exp(-1.5) rtol=0.1

        _, err = run_main_test("run walls.json dwi.json -o no_diff.csv --diffusivity 0. --bvecs=bvecs")
        no_diff = DataFrame(CSV.File("no_diff.csv"))
        @test size(no_diff, 1) == 2
        @test no_diff[!, :sequence] == [1, 1]
        @test no_diff[!, :bvec] == [1, 2]
        @test no_diff[1, :transverse] ≈ no_diff[1, :nspins]
        @test no_diff[2, :transverse] ≈ no_diff[2, :nspins]
    end
end
end
