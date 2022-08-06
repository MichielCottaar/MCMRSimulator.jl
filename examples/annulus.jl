import MRSimulator as mr
# Set up infinitely repeating aligned cylinders
outer = mr.Cylinder(0.9)
inner = mr.Cylinder(0.8)
geometry = mr.Repeated([inner, outer], [2., 2., Inf])

# Build spin echo diffusion-weighted sequence with perfect pulses and instantaneous gradients
TR = 2000.
TE = 80.
qval = 0.
sequence = mr.perfect_dwi(TR=TR, TE=TE, qval=qval)

micro = mr.Microstructure(diffusivity=mr.field(3.))

spins = [mr.Spin(position=rand(mr.PosVector)) for _ in 1:30]

@profview mr.evolve_TR(spins, sequence, micro)


res = mr.evolve_TR(spins, sequence, micro)