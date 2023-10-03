# [Tutorial](@id tutorial_cli)
This tutorial will walk through an example of modelling the MRI signal evolution for a diffusion-weighted sequence.
The spins in this simulation will be constrained by regularly packed cylinders.
This tutorial will use the command line interface, which we assume is available through the `mcmr` command (see [installation instructions](@ref installation)).
If you would prefer to interact with MCMRSimulator in Julia, you can find a tutorial doing the same analysis [here](@ref tutotial_julia).

In general, running a simulation will consist of the following three steps:
- Creating a geometry using one or more calls to `mcmr geometry` ([full description](@ref geometry))
- Defining a sequence using `mcmr sequence` ([full description](@ref sequence))
- Running the actual simulation using `mcmr run`
We will look through each of these steps below.

## Defining the geometry
First we will define a geometry formed of regularly packed axons.
This is represented by a single cylinder pointing in the z-direction with a radius of 1 micrometer that repeats itself every 2.5 micrometer (in both the x-, and y-direction).
```bash
mcmr geometry create cylinders 1 geometry.json --radius 1 --repeats 2.5 2.5
```
This will create a JSON file with the full information on the geometry:
```json
  {
     "type": "Cylinders",
     "number": 1,
    "#radius_description": "Radius of the cylinder. Field is required. Expected type: Float64.",
    "radius": 1.0,
    "#rotation_description": "Rotation applied to all obstructions in group. Can be set to a matrix or one of :x, :y, or, :z (see [`get_rotation`](@ref)). Field is required. Expected type: StaticArraysCore.SMatrix{3, 2, Float64, 6}.",
    "rotation": [[1.0,0.0,0.0],[0.0,1.0,0.0]],
    "#grid_resolution_description": "Resolution of the grid that the volume is split up into (um). Field is required. Expected type: Float64.",
    "grid_resolution": 1.0,
    "#R1_surface_description": "Additional longitudinal relaxation rate (kHz). Surface property. Field is required. Expected type: Float64.",
    "R1_surface": 0.0,
    "#R1_inside_description": "Additional longitudinal relaxation rate (kHz). Inside property. Field is required. Expected type: Float64.",
    "R1_inside": 0.0,
    "#R2_surface_description": "Additional transverse relaxation rate (kHz). Surface property. Field is required. Expected type: Float64.",
    "R2_surface": 0.0,
    "#R2_inside_description": "Additional transverse relaxation rate (kHz). Inside property. Field is required. Expected type: Float64.",
    "R2_inside": 0.0,
    "#off_resonance_surface_description": "Additional off-resonance field offset (kHz). Surface property. Field is required. Expected type: Float64.",
    "off_resonance_surface": 0.0,
    "#off_resonance_inside_description": "Additional off-resonance field offset (kHz). Inside property. Field is required. Expected type: Float64.",
    "off_resonance_inside": 0.0,
    "#position_description": "Spatial offset of obstruction from origin. Field is required. Expected type: StaticArraysCore.MVector{2, Float64}.",
    "position": [0.0,0.0],
    "#g_ratio_description": "Inner/outer radius used for susceptibility calculation Field can be null. Expected type: Float64.",
    "g_ratio": 1.0,
    "#susceptibility_iso_description": "Isotropic component of the susceptibility (in ppm). Field can be null. Expected type: Float64.",
    "susceptibility_iso": -0.1,
    "#susceptibility_aniso_description": "Ansotropic component of the susceptibility (in ppm). Field can be null. Expected type: Float64.",
    "susceptibility_aniso": -0.1,
    "#lorentz_radius_description": "Only compute field explicitly for cylinders with this Lorentz radius. Field can be null. Expected type: Float64.",
    "lorentz_radius": 5.0,
    "#repeats_description": "Length scale on which the obstructions are repeated (um). Field can be null. Expected type: StaticArraysCore.MVector{2, Float64}.",
    "repeats": [2.5,2.5],
    "#use_boundingbox_description": "Use bounding boxes for an initial filtering of possible intersections. Field can be null. Expected type: Bool.",
    "use_boundingbox": true,
    "#dwell_time_description": "Average time a particle stays stuck to the surface (ms). Surface property. Field can be null. Expected type: Float64.",
    "dwell_time": null,
    "#density_description": "Surface density of stuck particles relative to the volume density (um). Surface property. Field can be null. Expected type: Float64.",
    "density": null,
    "#permeability_description": "Probability of particle passing through the obstruction. Surface property. Field can be null. Expected type: Float64.",
    "permeability": null,
    "#relaxivity_description": "Fraction of transverse spin lost each time it hits the obstruction. Surface property. Field can be null. Expected type: Float64.",
    "relaxivity": null
  }
```
You can see how that the `repeats` and `radius` keywords have been set to our predefined values.
You can alter these and other geometry properties by editing this JSON directly or using the flags when creating the geometry.
For a full overview of these flags, you can run:
```
mcmr geometry create cylinders --help
```
How these various properties affect the simulation is described [here](@ref properties).

The procedure to create [`Walls`](@ref), [`Spheres`](@ref), or [`Annuli`](@ref) is very similar as for the [`Cylinders`](@ref) illustrated above.
Randomly distributed cylinders, annuli, and spheres can be created using `mcmr geometry create-random`.

## Defining the sequence
For now, the command line interface only supports
The next step is to define a sequence (see [here](@ref sequence) for more details). 
There are several built-in sequences available, which you can see listed by running:
```bash
mcmr sequence
```

Here, we will create a diffusion-weighted sequence:
```bash
mcmr sequence dwi dwi.json --bval=2 --TR=1000 --TE=80 --B0=3
```

This produces another JSON:
```json
{"scanner":{"B0":3.0,"gradient":null,"slew_rate":null},"gradients":[{"shape":{"times":[0.0,5.0e-324,39.99999999999999,40.0],"amplitudes":[[0.0,0.0,0.0],[0.0010896594058735262,0.0,0.0],[0.0010896594058735262,0.0,0.0],[0.0,0.0,0.0]]},"origin":[0.0,0.0,0.0]},{"shape":{"times":[40.0,40.00000000000001,79.99999999999999,80.0],"amplitudes":[[0.0,0.0,0.0],[0.0010896594058735262,0.0,0.0],[0.0010896594058735262,0.0,0.0],[0.0,0.0,0.0]]},"origin":[0.0,0.0,0.0]}],"instants":[{"time":0.0,"flip_angle":90.0,"cf":6.123233995736766e-17,"sf":1.0,"phase":-90.0,"cp":6.123233995736766e-17,"sp":-1.0},{"time":40.0,"flip_angle":180.0,"cf":-1.0,"sf":1.2246467991473532e-16,"phase":0.0,"cp":1.0,"sp":0.0}],"pulses":[],"TR":1000.0,"readout_times":[80.0]}
```
This one is less readable or editable by users, but basically describes the sequence diagram.

## Running the simulation
To get instructions on running the simulations, we can check the help message of `mcmr run`:
```bash
mcmr run --help
```

We can see that in addition to defining the geometry and the sequence, we can also control the simulation properties such as the `--diffusivity`, `--R1`, and `--R2`.

The simulation is initialised by randomly distributing a number of spins (set by `--Nspins`) uniformly across a bounding box with size given by `--voxel-size`.
This initial state might also contain bound spins (if the `--density` flag was set to a non-zero value during the geometry generation).

The DWI sequence defined above contains a [`Readout`](@ref) object at the echo time (80 ms). By default, this is used for readout:
```bash
mcmr run geometry.json dwi.json -o signal.csv
```
This produces the CSV file, which looks like
```csv
sequence,TR,readout,subset,nspins,longitudinal,transverse,phase,Sx,Sy
1,1,1,1,10000,-1.1262102361797588e-12,9997.399612964802,-0.006154308518773632,9997.399555292097,-1.0738501510661669
```

The columns in this file store the following information:
- "sequence": integer; index of the sequence (always 1 if only single sequence used)
- "TR": integer; index of the repetition time that this data was acquired (between 1 and the value of `--nTR`)
- "readout": integer; index of the readout within a TR.
- "subset": integer; index of the subset of the total signal (e.g., intra-axonal) that has been output (see the `--subset` flag). The total signal will always be included with "subset of 0.
- "time": time of the readout in ms (since start of simulation)
- "nspins": total number of spins contributing to the signal (might change for certain subsets of spins)
- "longitudinal": average longitudinal signal
- "transverse": average transverse signal
- "phase": average phase of the signal (in degrees)
- "Sx": signal strength in the x-direction
- "Sy": signal strength in the y-direction

A more complete state of all the spins can be produced using the `--output-snapshot` flag.
For example, the command
```bash
mcmr run geometry.json dwi.json -output-snapshot snapshot.csv
```
will produce a file named "snapshot.csv" with:
```csv
sequence,TR,readout,spin,x,y,z,longitudinal,transverse,phase,Sx,Sy
1,1,1,1,0.5049228156370124,0.05412808852066652,6.442542415589706,-2.220446049250313e-16,0.9999999999999942,0.25379615194714045,0.9999901894332786,0.004429563994801518
1,1,1,2,-0.9688485879911409,0.09305926325698886,-1.724633930537917,-2.220446049250313e-16,0.9999999999999917,-0.0047769059448796725,0.999999996524485,-8.337273670013147e-5
1,1,1,3,0.2629314013825268,-0.6671294575623549,-17.506721187691674,-2.220446049250313e-16,0.9999999999999946,0.6082085693901605,0.9999436588469045,0.010615042715627557
1,1,1,4,-0.9924343175166189,-0.051930346020524064,-11.845895413220594,0.0,0.9999999999999928,-0.11879290332484516,0.9999978506577348,-0.0020733258055911637
1,1,1,5,0.08578653895854582,0.5766674996691633,2.905211448847795,-2.220446049250313e-16,0.9999999999999931,-0.28110585975525737,0.9999879645130568,-0.004906203116436347
1,1,1,6,0.5283705886599059,-0.6771241777295106,18.97063962214857,-2.220446049250313e-16,0.9999999999999921,-1.2036996551270818,0.999779329114167,-0.021006976841285743
1,1,1,7,0.2676027006707101,0.27872400783545637,21.212386477239814,0.0,0.999999999999992,1.8741387410642074,0.9994650791902853,0.03270405906215431
...
```
Each row corresponds to the state of a single spin. In addition to all the columns listed above, we now have 4 more columns:
- "spin": integer; index of the spin
- "x"/"y"/"z": floats; position of the spin at the time of the readout

The readout times can be adjusted using the `--nTR`, `--time`, and `--skip-TR` flags.
For more examples of this, see the end of the [tutorial using the Julia interface](@ref tutorial_julia).
