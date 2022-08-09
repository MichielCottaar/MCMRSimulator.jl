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

snap = mr.Snapshot(3000);

simulation = mr.Simulation(snap, [sequence], micro);

append!(simulation, 2.);
@time append!(simulation, 200.);
@profview append!(simulation, 200.);

before = mr.isinside(simulation.latest[1], geometry);
after = mr.isinside(simulation.latest[end], geometry);
t = simulation.latest[end].time
frac = sum(xor.(before, after)) / length(snap.spins)
println("Fraction of spins crossing boundary in $t ms: $frac")


function test_dwi(;nspins=100000, bval=2., store_every=1., diffusivity=3.)
    snap = mr.Snapshot(nspins)
    sequence = mr.perfect_dwi(bval=bval)
    sim = mr.Simulation(snap, [sequence], diffusivity=diffusivity, store_every=store_every)
    append!(sim, 90)
    span = sim.readout[1][1]
    log(mr.transverse(span) / nspins)
end