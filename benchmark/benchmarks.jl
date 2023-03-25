using BenchmarkTools
import MCMRSimulator as mr


span = mr.Snapshot(300)
sequence = mr.dwi(bval=2.)
geometries = (
    mr.spheres(0.9, repeats=[2, 2, 2]),
    mr.cylinders(0.9, repeats=[2, 2]),
    mr.walls(repeats=2),
    mr.TransformObstruction(mr.box_mesh(), repeats=[1.5, 1.5, 1.5]),
)


SUITE = BenchmarkGroup()
SUITE["no diffusion"] = @benchmarkable mr.evolve(span, mr.Simulation([sequence], diffusivity=0.), sequence.TR)
SUITE["Repeating spheres"] = @benchmarkable mr.evolve(span, mr.Simulation([sequence], diffusivity=3., geometry=$geometries[1]), sequence.TR)
SUITE["Repeating cylinders"] = @benchmarkable mr.evolve(span, mr.Simulation([sequence], diffusivity=3., geometry=$geometries[2]), sequence.TR)
SUITE["Repeating walls"] = @benchmarkable mr.evolve(span, mr.Simulation([sequence], diffusivity=3., geometry=$geometries[3]), sequence.TR)
SUITE["Repeating mesh boxes"] = @benchmarkable mr.evolve(span, mr.Simulation([sequence], diffusivity=3., geometry=$geometries[4]), sequence.TR)