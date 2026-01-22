"""
Oceananigans based model of an ocean gyre in a bounded domain.

## Details

Model of wind and buoyancy forced ocean gyre adapted from MITgcm documentation
[baroclinic ocean gyre
example](https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html)
implemented using Oceananigans (Ramadhan et al. 2020).

Compared to MITgcm example, beyond the change in the underlying software framework, some
key differences in the (default) model configuration are

- The model uses a seawater buoyancy formulation with tracer fields corresponding to
  salinity and (potential) temperature used to determine the buoyancy using a non-linear
  (polynomial) equation of state (Roquet et al. 2015).
- A weighted essentially non-oscillatory numerical scheme is used for momentum advection
  (Silvestri et al. 2024).
- The turbulence closure / parameterization used for small-scale turbulent mixing is a
  convective adjustment turbulent kinetic energy (CATKE) scheme (Wagner et al. 2025).
- Adaptive time stepping is used based on controlling advective CFL number.

$(EXPORTS)

## References

1. MITgcm contributors (2025). MITgcm/MITgcm: checkpoint69e (Version checkpoint69e).
   Zenodo. https://doi.org/10.5281/zenodo.15320163
2. Ramadhan, A., Wagner, G., Hill, C., Campin, J. M., Churavy, V., Besard, T.,
   Souza, A., Edelman, A., Ferrari, R. & Marshall, J. (2020).
   Oceananigans. jl: Fast and friendly geophysical fluid dynamics on GPUs.
   Journal of Open Source Software, 5(53).
3. Roquet, F.; Madec, G.; Brodeau, L. and Nycander, J. (2015).
   Defining a simplified yet “realistic” equation of state for seawater.
   Journal of Physical Oceanography 45, 2564-2579.
4. Silvestri, S., Wagner, G. L., Campin, J. M., Constantinou, N. C., Hill, C. N., 
   Souza, A., & Ferrari, R. (2024).
   A new WENO-based momentum advection scheme for simulations of ocean mesoscale
   turbulence. Journal of Advances in Modeling Earth Systems, 16(7), e2023MS004130.
5. Wagner, G. L., Hillier, A., Constantinou, N. C., Silvestri, S., Souza, A.,
   Burns, K. J., Hill, C., Campin, J-M., Marshall, J. & Ferrari, R. (2025).
   Formulation and calibration of CATKE, a one-equation parameterization for microscale
   ocean mixing. Journal of Advances in Modeling Earth Systems, 17(4), e2024MS004522.
"""
module GyreInABox

using DocStringExtensions
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState
using CairoMakie
using Printf

export GyreInABoxParameters, GyreInABoxConfiguration
export HorizontalSlice, LongitudeDepthSlice, LatitudeDepthSlice
export DepthAveraged, FreeSurfaceFields, MOCStreamFunction, BarotropicStreamFunction
export setup_model, initialize_model!, setup_simulation
export run_simulation, record_animations

include("parameters.jl")
include("grids.jl")
include("boundary_conditions.jl")
include("initial_conditions.jl")
include("models.jl")
include("simulations.jl")
include("outputs.jl")

end
