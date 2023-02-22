"""
    Scanner(;B0=3., max_gradient=Inf, max_slew_rate=Inf, units=:kHz)

Properties of an MRI scanner relevant for the MR signal simulations.
- [`B0`](@ref): magnetic field strength (in Tesla)
- [`max_gradient`](@ref): maximum gradient strength long each axis.
- [`max_slew_rate`](@ref): maximum rate of change in the gradient strength

By default `gradient` and `slew_rate` are expected to be provided in units of, respectively, kHz/um and kHz/um/ms.
However, if the keyword `units=:Tesla` is set, the `gradient` and `slew_rate` should be provided in units of, respectively, mT/m and T/m/s.
"""
struct Scanner
    B0::Float
    gradient::Float
    slew_rate::Float
    function Scanner(;B0=3., gradient=Inf, slew_rate=Inf, units=:kHz)
        if isinf(gradient) && !isinf(slew_rate)
            error("Can't have infinite gradient strength with finite slew rate.")
        end
        if units == :Tesla
            gradient *= gyromagnetic_ratio * 1e-9
            slew_rate *= gyromagnetic_ratio * 1e-9
        end
        new(Float(B0), Float(gradient), Float(slew_rate))
    end
end

"""
    B0(scanner)
    B0(sequence)

Returns the magnetic field strength of the scanner in Tesla.
"""
B0(scanner::Scanner) = scanner.B0

"""
    max_gradient(scanner[, units])

Returns the maximum magnetic field gradient of the scanner in kHz/um.
By setting `units` to :Tesla, the gradient strength can be returned in mT/m instead.
"""
max_gradient(scanner::Scanner, units=:kHz) = units == :kHz ? scanner.gradient : scanner.gradient / (gyromagnetic_ratio * 1e-9)

"""
    max_slew_rate(scanner[, units])

Returns the maximum magnetic field slew rate of the scanner in kHz/um/ms.
By setting `units` to :Tesla, the slew rate can be returned in T/m/s instead.
"""
max_slew_rate(scanner::Scanner, units=:kHz) = units == :kHz ? scanner.slew_rate : scanner.slew_rate / (gyromagnetic_ratio * 1e-9)


"""
Siemens MAGNETOM 3T Prisma MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/3t-mri-scanner/magnetom-prisma).
"""
Siemens_Prisma = Scanner(B0=3., gradient=80, slew_rate=200, units=:Tesla)

"""
Siemens MAGNETOM 7T Terra MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/7t-mri-scanner/magnetom-terra)
"""
Siemens_Terra = Scanner(B0=7., gradient=80, slew_rate=200, units=:Tesla)

"""
Siemens 3T Connectom MRI scanner ([fan22_MappingHumanConnectome](@cite)).
"""
Siemens_Connectom = Scanner(B0=3., gradient=300, slew_rate=200, units=:Tesla)