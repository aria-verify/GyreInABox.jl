```@meta
CurrentModule = GyreInABox
```

# GyreInABox

Documentation for [GyreInABox](https://github.com/aria-verify/GyreInABox.jl).

## Model details

The model is an adaptation of the [baroclinic gyre
example](https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html)
from the MITgcm documentation. It simulates a wind and buoyancy forced double-gyre ocean
circulation on a spherical shell sector spatial domain, with default domain extents
$(0^\circ, 60^\circ) \times (15^\circ, 75^\circ) \times (-1800\textsf{m}, 0\textsf{m})$
(longitude × latitude × depth dimensions).

By default a non-linear [TEOS-10 polynomial equation of
state](https://clima.github.io/SeawaterPolynomials.jl/stable/#The-TEOS-10-standard) is
used to compute buoyancy from salinity and temperature fields. No slip / no-flux
boundary conditions are applied to the velocity fields on all walls, and a boundary
condition corresponding to a damping zonal drag on the bottom surface. The zonal
velocity component is subject to a latitude dependent wind stress roughly reflecting
average zonal wind patterns at the surface. The temperature and salinity fields have
no-flux boundary conditions applied on all walls and the bottom surface, and relaxation
boundary conditions on the top surface which restore the surface fields towards a
latitude dependent reference temperature / salinity with a parameterized relaxation
time. The surface wind stress, reference temperature and reference salinity as a
function of latitude are shown in the figure below.

```@setup model_details
using GyreInABox
using CairoMakie
using Oceananigans

parameters = GyreInABoxParameters()
configuration = GyreInABoxConfiguration()

grid = GyreInABox.setup_grid(configuration)

λ, φ, z = nodes(grid, Center(), Center(), Center())

wind_stress = GyreInABox.zonal_wind_stress.(Nothing, φ, Nothing, (parameters,))
reference_temperature =  GyreInABox.reference_surface_temperature.(φ, (parameters,))
reference_salinity =  GyreInABox.reference_surface_salinity.(φ, (parameters,))

figure = Figure(size=(800, 200))

ylabel = "Latitude φ / °"

axis_wind_stress = Axis(figure[1, 1]; ylabel, xlabel="Zonal wind stress / 10⁻³ m² s⁻²")
axis_temperature = Axis(figure[1, 2]; ylabel, xlabel="Reference temperature / °C")
axis_salinity = Axis(figure[1, 3]; ylabel, xlabel="Reference salinity /  g kg⁻¹")

# Due to Oceananigans' flux convention, positive zonal wind stress flux at top surface
# boundary produces currents in negative (zonal) direction and vice versa therefore
# negate wind stress profile to better match visual intuition
lines!(axis_wind_stress, -wind_stress * 10^3, φ)
lines!(axis_temperature, reference_temperature, φ)
lines!(axis_salinity, reference_salinity, φ)

save("surface_forcings.svg", figure)
```

![Surface forcings for model](surface_forcings.svg)

The model is initialised from rest (zero velocity). The temperature and salinity fields
are initialised to match the reference temperature and salinity at the surface and
following smoothly varying thermocline / halocline like profiles with depth, and
constant in longitude, as shown in the figure below.

```@setup model_details
initial_temperature = GyreInABox.initial_temperature.(φ, z', (parameters,))
initial_salinity = GyreInABox.initial_salinity.(φ, z', (parameters,))

figure = Figure(size=(800, 300))

xlabel = "Latitude φ / °"
ylabel = "Depth / m"

axis_temperature = Axis(figure[1, 1]; xlabel, ylabel, title="Initial temperature")
axis_salinity = Axis(figure[1, 3]; xlabel, ylabel, title="Initial salinity")

hm_T = heatmap!(axis_temperature, φ, z, initial_temperature; colormap=:thermal)
Colorbar(figure[1, 2], hm_T; label="ᵒC")
hm_S = heatmap!(axis_salinity, φ, z, initial_salinity; colormap=:haline)
Colorbar(figure[1, 4], hm_S; label="g / kg")

save("initial_temperature_and_salinity.svg", figure)
```

![Initial temperature and salinity fields for model](initial_temperature_and_salinity.svg)

## Usage example

The below example illustrates setting up and running a simulation of the model with the
default parameter values and configuration modulo changing the model to record only a
horizontal (surface) slice and depth slice along longitude axis as outputs. The
resulting horizontal and depth slices of simulated fields are recorded as an animation
using CairoMakie. A coarse default grid size (60 × 60 × 15) and small simulation time
(60 days) are set so as to allow a simulation to run in around 5 minutes on a CPU, but
if running on a GPU (controlled by `architecture` keyword argument to
`GyreInABoxConfiguration`) then finer spatial discretizations and longer simulation
times can easily be used.

```@example
using Oceananigans.Units
using GyreInABox

parameters = GyreInABoxParameters()
output_types = (LongitudeDepthSlice(), HorizontalSlice())
configuration = GyreInABoxConfiguration(output_types=output_types)
run_simulation(parameters, configuration)
for output_type in output_types
    record_animation(configuration, output_type)
end
nothing # hide
```

![Example simulated output (longitude-depth slice)](gyre_model_longitude_depth_slice.mp4)
![Example simulated output (horizontal slice)](gyre_model_horizontal_slice.mp4)

## API reference

```@index
```

```@autodocs
Modules = [GyreInABox]
```
