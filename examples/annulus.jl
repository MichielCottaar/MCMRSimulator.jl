import MRSimulator as mr
using StaticArrays
# Set up infinitely repeating aligned cylinders
outer = mr.Cylinder(0.9)
inner = mr.Cylinder(0.8)
geometry = mr.Repeated([inner, outer], [2., 2., Inf])

# Build spin echo diffusion-weighted sequence with perfect pulses and instantaneous gradients
TR = 200.
TE = 80.
bval = 2.
sequence = mr.perfect_dwi(TR=TR, TE=TE, bval=bval)

micro = mr.Microstructure(diffusivity=mr.field(3.), geometry=geometry)

spins = [mr.Spin(position=SA_F64[0., 0., z]) for z in 1:30000];

simulation = mr.Simulation(spins, [sequence], micro);

append!(simulation, 2.);
@time append!(simulation, 20.);
@profview append!(simulation, 200.);
