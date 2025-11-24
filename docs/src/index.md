```@meta
CurrentModule = GyreInABox
```

# GyreInABox

Documentation for [GyreInABox](https://github.com/aria-verify/GyreInABox.jl).

## Usage example

The below example illustrates setting up and running a simulation of the model with the
default parameter values, and default configuration modulo reduction of spatial
resolution to $60\times 60 \times 25$ and the simulation time to 6 hours and changing
the model to record only a depth slice along longitude axis as outputs. The resulting
depth slices of simulated fields are recorded as an animation using CairoMakie.

```@example

using Oceananigans.Units
using GyreInABox

parameters = GyreInABoxParameters()
output_type = LongitudeDepthSlice()
configuration = GyreInABoxConfiguration(
    grid_size=(60, 60, 25),
    simulation_time=6hour,
    output_types=(output_type,)
)
run_simulation(parameters, configuration)
record_animation(configuration, output_type)
nothing # hide
```

![Example simulated output](gyre_model_longitude_depth_slice.mp4)

## API reference

```@index
```

```@autodocs
Modules = [GyreInABox]
```
