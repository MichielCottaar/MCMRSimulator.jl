using MCMRSimulator
using MRIBuilder
pos, rad = random_positions_radii(30., 0.6, 2; variance=0.)
main_seq = DWI(bval=3., slice_thickness=2.)
geometries = Dict(
    "cylinders" => Cylinders(position=pos, radius=rad, repeats=(30., 30.)),
    "cylinders_myelin" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), g_ratio=0.8),
    "cylinders_MT" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), surface_density=1., dwell_time=30.),
    "cylinders_perm" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), permeability=1.),
)

simulations = Dict{String, Simulation}(
    key => Simulation(main_seq; geometry=geometry, diffusivity=2.)
    for (key, geometry) in pairs(geometries)
)
simulations["no_diff"] = Simulation(main_seq; diffusivity=0.);