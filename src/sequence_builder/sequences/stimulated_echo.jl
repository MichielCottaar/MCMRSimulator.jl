"""Stimulated Echo implementation"""
module stimulated_echo
import ....Scanners: Scanner
import ....Sequences: InstantRFPulse, Readout, RFPulse, InstantGradient
import ...DefineSequence: define_sequence
import ...Diffusion: add_linear_diffusion_weighting
import ...BuildingBlocks: duration
import StaticArrays: SVector

end