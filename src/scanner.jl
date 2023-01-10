"""
    Scanner(;B0=3., gradient=Inf, slew_rate=Inf)

Properties of an MRI scanner relevant for the MR signal simulations.
- `B0`: magnetic field strength (in Tesla)
- `gradient`: maximum gradient strength (in mT/m) long each axis
- `slew_rate`: maximum rate of change in the gradient strength (in T/m/s)
"""
struct Scanner
    B0::Float
    gradient::Float
    slew_rate::Float
    function Scanner(;B0=3., gradient=Inf, slew_rate=Inf)
        if isinf(gradient) && !isinf(slew_rate)
            error("Can't have infinite gradient strength with finite slew rate.")
        end
        new(Float(B0), Float(gradient), Float(slew_rate))
    end
end


"""
Siemens MAGNETOM 3T Prisma MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/3t-mri-scanner/magnetom-prisma).
"""
Siemens_Prisma = Scanner(B0=3., gradient=80, slew_rate=200)

"""
Siemens MAGNETOM 7T Terra MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/7t-mri-scanner/magnetom-terra)
"""
Siemens_Terra = Scanner(B0=7., gradient=80, slew_rate=200)

"""
Siemens 3T Connectom MRI scanner ([fan22_MappingHumanConnectome](@cite)).
"""
Siemens_Connectom = Scanner(B0=3., gradient=300, slew_rate=200)