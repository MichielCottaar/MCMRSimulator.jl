import MRSimulator as mr
using StaticArrays
# Set up infinitely repeating aligned cylinders
outer = mr.Cylinder(0.9)
inner = mr.Cylinder(0.8)
geometry = mr.Repeated([inner, outer], [2., 2., Inf])

mb = mr.Repeated(mr.box_mesh(grid_size=3), [1.5, 1.5, 1.5])

# Build spin echo diffusion-weighted sequence with perfect pulses and instantaneous gradients
TR = 200.
TE = 80.
bval = 2.
sequence = mr.perfect_dwi(TR=TR, TE=TE, bval=bval)

snap = mr.Snapshot(3000);

simulation = mr.Simulation(snap, [sequence], diffusivity=3., geometry=[mb]);

append!(simulation, 2.);
@time append!(simulation, 200.);
@profview append!(simulation, 200.);