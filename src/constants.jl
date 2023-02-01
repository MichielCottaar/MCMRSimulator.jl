"gyromagnetic ratio of a proton (1H) in water in kHz/T"
#const gyromagnetic_ratio = 0.2675223151  # (10^6 rad⋅ms^−1⋅T^−1)
const gyromagnetic_ratio = 42576.38476  # (kHz/T)

"Float type used during the simulation (default: Float64). This has to be changed in the code."
const Float = Float64
const SA = SA_F64

"General definition used for length-3 vectors. Mostly used for positions."
const PosVector = SVector{3, Float}