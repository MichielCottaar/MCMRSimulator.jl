import MRSimulator as mr
using StaticArrays
# Set up infinitely repeating aligned cylinders
geometry = mr.cylinders([0.8, 0.9], shifts=[[0, 0], [0, 0]], repeats=[2., 2.])

#geometry = mr.TransformObstruction(mr.box_mesh(grid_size=10), repeats=[1.5, 1.5, 1.5])

# Build spin echo diffusion-weighted sequence with perfect pulses and instantaneous gradients
TR = 200.
TE = 80.
bval = 2.
sequence = mr.perfect_dwi(TR=TR, TE=TE, bval=bval)

snap = mr.Snapshot(3000);

simulation = mr.Simulation(snap, [sequence], diffusivity=3., geometry=geometry);

append!(simulation, 2.);
@time append!(simulation, 20.);
@profview append!(simulation, 20.);