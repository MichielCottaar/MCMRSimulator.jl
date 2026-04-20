# This script required MRIBuilder v0.3.0
# It will only work with Julia v1.11, not later versions
# This script generates the .seq files used in the tests.
# By pregenerating these .seq files, we can run the tests without MRIBuilder installed.

using MRIBuilder

for TE in [2, 20, 30, 100, 1000, 1e5]
    sequence = GradientEcho(TE=TE)
    write_sequence(joinpath(@__DIR__, "gradient_echo_TE_$(Int(TE)).seq"), sequence)
end

sequence_top = GradientEcho(TE=2.8)
write_sequence(joinpath(@__DIR__, "gradient_echo_TE_2.8.seq"), sequence_top)

for (sequence, filename) in [
    (DWI(TE=80., bval=0.5), "dwi_te_80_bval_0.5.seq"),
    (DWI(TE=80., bval=1), "dwi_te_80_bval_1.seq"),
    (DWI(TE=80., bval=2), "dwi_te_80_bval_2.seq"),
    (DWI(bval=0.1, TE=80.02, Δ=40.03, gradient=(type=:instant, )), "dwi_te_80_bval_0.1_instant_diffusion_time_40.seq"),
    (DWI(TE=80, bval=2., gradient=(type=:instant, ), diffusion_time=40.), "dwi_te_80_bval_2_instant_diffusion_time_40.seq"),
    (DWI(TE=80, bval=2., gradient=(type=:instant, )), "dwi_te_80_bval_2_instant.seq"),
    (DWI(bval=2, TE=80, scanner=Siemens_Prisma, gradient=(rise_time=:min,)), "dwi_te_80_bval_2_prisma_min_rise_time.seq"),
]
    write_sequence(joinpath(@__DIR__, filename), sequence)
end

for (δ, Δ) in [
    (nothing, 30.),
    (10, 30.),
    (10., 30.),
    (20, 40.),
    (nothing, nothing),
    (10., nothing),
    (10, nothing),
    (0, nothing),
    (0, 70.),
    (0, 0.1),
    (0, 70),
    (0, 30.),
]
    gradient = (isnothing(δ) || !iszero(δ)) ? (type=:trapezoid, δ=δ) : (type=:instant, )
    sequence = DWI(bval=0.3, TE=80, Δ=Δ, gradient=gradient, scanner=Siemens_Connectom)
    write_sequence(joinpath(@__DIR__, "dwi_te_80_bval_0.3_gradient_δ_$(δ)_Δ_$(Δ).seq"), sequence)
end

for qval in [0.01, 0.1, 1.]
    sequence = DWI(TE=101, diffusion_time=100, gradient=(type=:instant, qval=qval))
    write_sequence(joinpath(@__DIR__, "dwi_te_101_diffusion_time_100_qval_$(qval).seq"), sequence)

    sequence2 = DWI(TE=101, diffusion_time=100, gradient=(type=:instant, orientation=[1., 0., 0.], qval=qval))
    write_sequence(joinpath(@__DIR__, "dwi_te_101_diffusion_time_100_qval_$(qval)_orientation_x.seq"), sequence2)
end

for dt in [0.003, 0.01]
    sequence = DWI(diffusion_time=dt, TE=2, bval=2., gradient=(type=:instant, ))
    write_sequence(joinpath(@__DIR__, "dwi_te_2_bval_2_diffusion_time_$(dt).seq"), sequence)
end