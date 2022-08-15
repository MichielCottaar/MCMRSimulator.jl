using BenchmarkTools
import MRSimulator as mr


span = mr.Snapshot(300)
sequence = mr.perfect_dwi(bval=2.)
geometries = (
    [mr.Repeated(mr.Sphere(0.9), [2., 2., 2.])],
    [mr.Repeated(mr.Cylinder(0.9), [2., 2., Inf])],
    [mr.Repeated(mr.Wall(), [2., Inf, Inf])],
    [mr.Repeated(mr.box_mesh(), [1.5, 1.5, 1.5])],
)


SUITE = BenchmarkGroup()
SUITE["no diffusion"] = @benchmarkable append!(mr.Simulation(span, [sequence], diffusivity=0.), sequence.TR)
SUITE["Repeating spheres"] = @benchmarkable append!(mr.Simulation(span, [sequence], diffusivity=3., geometry=$geometries[1]), sequence.TR)
SUITE["Repeating cylinders"] = @benchmarkable append!(mr.Simulation(span, [sequence], diffusivity=3., geometry=$geometries[2]), sequence.TR)
SUITE["Repeating walls"] = @benchmarkable append!(mr.Simulation(span, [sequence], diffusivity=3., geometry=$geometries[3]), sequence.TR)
SUITE["Repeating mesh boxes"] = @benchmarkable append!(mr.Simulation(span, [sequence], diffusivity=3., geometry=$geometries[4]), sequence.TR)