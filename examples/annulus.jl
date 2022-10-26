import MRSimulator as mr
using StaticArrays
# Set up infinitely repeating aligned cylinders
@time (positions, radii) = mr.random_positions_radii((30., 30.), 0.7, 2);
geometry = mr.annuli(0.8 .* radii, radii, positions=positions, repeats=[30., 30.], MT_fraction=1e-3, orientation=:y, myelin=true)

#geometry = mr.TransformObstruction(mr.box_mesh(grid_size=10), repeats=[1.5, 1.5, 1.5])

# Build spin echo diffusion-weighted sequence with perfect pulses and instantaneous gradients
TR = 200.
TE = 80.
bval = 2.
sequence = mr.perfect_dwi(TR=TR, TE=TE, bval=bval)

snap = mr.Snapshot(3000);

simulation = mr.Simulation([sequence], diffusivity=3., geometry=geometry);

mr.evolve(snap, simulation, 2.);
@time mr.evolve(snap, simulation, 200.);
@profview mr.evolve(snap, simulation, 200.);