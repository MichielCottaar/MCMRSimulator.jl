@testset "test_cli.jl" begin
import MCMRSimulator.CLI: run_main_test
import JSON: parse
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
                _, err = run_main_test("geometry create cylinders 2 test.json --radius 1 2 --position 0 1 2 3 --R1_surface 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
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
                _, err = run_main_test("geometry create spheres 2 test.json --radius 1 2 --position 0 1 2 3 4 5 --R1_surface 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
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
                _, err = run_main_test("geometry create walls 2 test.json --position 0 1 --R1 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Walls"
                @test result["number"] == 2
                @test result["position"] == [0., 1.]
                @test result["R1"] == 1.5
                @test result["R2"] == 0.
            end
        end
        @testset "create annuli" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create annuli 2 test.json --inner 1 2 --outer 3 --myelin")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test result["number"] == 2
                @test result["position"] == [0., 0.]
                @test result["inner"] == [1, 2]
                @test result["outer"] == 3.
                @test result["myelin"] == true
            end
        end
    end
    @testset "mcmr geometry create-random" begin
        @testset "random cylinders" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random cylinders 0.5 test.json --repeats 20 20")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random spheres" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random spheres 0.5 test.json --repeats 20 20 20")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random annuli" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random annuli 0.5 test.json --repeats 20 20 --g-ratio=0.8")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
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
                @test length(sequence.instants) == 1
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
                @test length(sequence.instants) == 2
                @test sequence.instants[1].time == 0.
                @test sequence.instants[2].time == 15.
                @test length(sequence.gradients) == 0
                @test sequence.TR == 100
                @test sequence.readout_times == [30.]
            end
        end
        @testset "spin_echo with finite pulses" begin
            in_tmpdir() do
                _, err = run_main_test("sequence spin-echo seq.json --TE 30 --TR 100 --excitation-duration 2 --refocus-duration 4")
                @test length(err) == 0
                sequence = mr.read_sequence("seq.json")
                @test length(sequence.pulses) == 2
                @test length(sequence.instants) == 0
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
                @test length(sequence.instants) == 2
                @test sequence.instants[1].time == 0.
                @test sequence.instants[2].time == 15.
                @test length(sequence.gradients) == 2
                @test sequence.TR == 100
                @test sequence.readout_times == [30.]
            end
        end
    end
    @testset "Test full simulations" begin
        @testset "Setting R2" begin
            in_tmpdir() do
                @testset "Single global R2" begin
                    _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2 2.2 2.2")
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
                    _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2 2.2 2.2 --R2_inside=0.2")
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
                    @test all(result[!, :sequence] .== 1)
                end
            end
        end
    end
    @testset "Setting R1" begin
        _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2 2.2 2.2 --R1_inside=0.02")
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
        @test all(result[!, :sequence] .== 1)
    end
end

end
