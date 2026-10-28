# Wind, temperature and salinity forced gyre with basin bathymetry

## Model details

This idealized model of the subpolar gyre has surface wind, temperature and salinity forcing and simplified topography with a northern basin separated from a southern open ocean region, based on the model described in [Spall (2011)](https://doi.org/10.1175/2011JCLI4130.1) and [Spall (2012)](https://doi.org/10.1175/JPO-D-11-0230.1).

```@setup model_details
using GyreInABox
using CairoMakie

parameters = Spall2011Parameters()
figure = GyreInABox.plot_domain_and_forcing(
    parameters, axis_height=450, axis_width=300, font_size=14
)
save("spall_2011_setup.svg", figure)
```

![Idealized subpolar gyre model bathymetry and surface forcing](spall_2011_setup.svg)

Spall (2011) describes an idealized marginal-sea configuration to represent the circulation which takes place in the subpolar North Atlantic Ocean and Nordic Seas. A northern basin is separated from a southern open-ocean region by a sill with bathymetry as shown in leftmost panel above, with domain spanning 1000 km zonally and 2000 km meridionally. The sea surface temperature is restored to a spatially varying profile, shown in second panel, loosely representing average surface atmospheric temperatures which are cooler in the north and northwest and warmer in the south. Surface wind forcing is applied with a meridionally varying wind stress profile corresponding to zonal velocity flux shown in the third panel panel. A linear equation of state is used with density linearly depending on temperature. The temperature in a sponge region at the south of the domain is restored to give a uniform vertical density gradient (stratification), approximately representing the thermodynamic processes maintaining stratification at lower latitudes. This restoration in the southern sponge region along with inclusion of a convective adjustment turbulence closure scheme enhancing vertical diffusivity for statically unstable conditions and shallow domain depth all contribute to a relatively short time to reach statistical equilibrium of order 50 simulated years from an at-rest initialization.

Spall (2012) adopts the same basic idealized configuration as Spall (2011) but with the addition of a separate salinity representation and freshwater forcing over the northern basin. The equation of state is modified to include a linear dependence on both temperature and salinity. A virtual surface salinity flux representing an influx of freshwater into the northern region is applied over the basin to the north of the sill, with a net evaporation-precipitation flux applied in this northern region, and zero elsewhere, as illustrated in the fourth pane above. In the southern sponge region the salinity is restored to a constant reference value with the same relaxation timescale as the temperature restoration.

## Usage example

The below example illustrates setting up and running a simulation of the model with the
default parameter values other than a reduced grid dimensionality. The resulting horizontal and depth slices of simulated fields are recorded as an animation using CairoMakie.

```@example
using Oceananigans.Units
using GyreInABox

parameters = Spall2011Parameters(grid_size=(50, 100, 25))
configuration = SimulationConfiguration(
    simulation_time=60days,
    output_filename_stem = "spall_2011_gyre_model",
    output_types=(
        y_depth_slice_output(; x_or_longitude=500kilometer),
        horizontal_slice_output(; depth=0.)
    ),
    progress_message_interval=1000
)
run_simulation(parameters, configuration)
plot_outputs(AnimationPlotOutput(), parameters, configuration)
nothing # hide
```

![Example simulated output (y-depth slice)](spall_2011_gyre_model_spatial_slice_at_x_500000.0_of_variables_u_v_w_S_T_at_time_interval_1day.mp4)
![Example simulated output (horizontal slice)](spall_2011_gyre_model_spatial_slice_at_z_0.0_of_variables_u_v_w_S_T_at_time_interval_1day.mp4)
