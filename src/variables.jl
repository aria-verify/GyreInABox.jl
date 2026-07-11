
abstract type AbstractModelVariable end

"""
    $(FUNCTIONNAME)(variable, model)

Field from `model` associated with `variable`.
"""
function field end

"""
    $(FUNCTIONNAME)(variable)

Short name associated with `variable` used as key in output files.
"""
function short_name end

"""
    $(FUNCTIONNAME)(variable)

CF conventions standard name associated with `variable`.
"""
function standard_name end

"""
    $(FUNCTIONNAME)(variable)

Human-readable long name associated with `variable` for use in for example plot labels.
"""
function long_name end

"""
    $(FUNCTIONNAME)(variable)

Physical units associated with `variable`.
"""
function units end

"""
    $(FUNCTIONNAME)(variable)

Color map from `cmocean` set to use in visualizations of `variable`.
"""
function color_map end

"""
    $(FUNCTIONNAME)(variable)

Spatial dimensions of field associated with `variable`.
"""
function spatial_dimensions(::AbstractModelVariable)
    SpatialDimensions(
        FreeSpatialDimension(), FreeSpatialDimension(), FreeSpatialDimension()
    )
end

"""
    $(FUNCTIONNAME)(variable)

Dictionary of any additional attributes to record in NetCDF output for `variable`.
"""
additional_attributes(::AbstractModelVariable) = Dict{String,String}()

"""
    $(FUNCTIONNAME)(variable)

Whether to symmetrize automatic value limits around zero for visualizations of `variable`.
"""
symmetrize_limits(::AbstractModelVariable) = true

abstract type AbstractVelocityVariable <: AbstractModelVariable end

units(::AbstractVelocityVariable) = "m s⁻¹"
color_map(::AbstractVelocityVariable) = :balance

abstract type AbstractSpecificKineticEnergyVariable <: AbstractModelVariable end

units(::AbstractSpecificKineticEnergyVariable) = "m² s⁻²"
color_map(::AbstractSpecificKineticEnergyVariable) = :amp
symmetrize_limits(::AbstractSpecificKineticEnergyVariable) = false

abstract type AbstractVelocityStreamFunctionVariable <: AbstractModelVariable end

units(::AbstractVelocityStreamFunctionVariable) = "m³ s⁻¹"
color_map(::AbstractVelocityStreamFunctionVariable) = :balance

abstract type AbstractBarotropicVelocityVariable <: AbstractModelVariable end

units(::AbstractBarotropicVelocityVariable) = "m² s⁻¹"
color_map(::AbstractBarotropicVelocityVariable) = :balance
function spatial_dimensions(::AbstractBarotropicVelocityVariable)
    SpatialDimensions(
        FreeSpatialDimension(), FreeSpatialDimension(), ReducedSpatialDimension()
    )
end

"""
$(TYPEDEF)

Eastward (zonal or x) component of velocity.
"""
struct EastwardVelocity <: AbstractVelocityVariable end

field(::EastwardVelocity, model) = model.velocities.u
short_name(::EastwardVelocity) = "u"
standard_name(::EastwardVelocity) = "eastward_sea_water_velocity"
long_name(::EastwardVelocity) = "Eastward velocity"

"""
$(TYPEDEF)

Northward (meridional or y) component of velocity.
"""
struct NorthwardVelocity <: AbstractVelocityVariable end

field(::NorthwardVelocity, model) = model.velocities.v
short_name(::NorthwardVelocity) = "v"
standard_name(::NorthwardVelocity) = "northward_sea_water_velocity"
long_name(::NorthwardVelocity) = "Northward velocity"

"""
$(TYPEDEF)

Upward (vertical or z) component of velocity.
"""
struct UpwardVelocity <: AbstractVelocityVariable end

field(::UpwardVelocity, model) = model.velocities.w
short_name(::UpwardVelocity) = "w"
standard_name(::UpwardVelocity) = "upward_sea_water_velocity"
long_name(::UpwardVelocity) = "Upward velocity"

"""
$(TYPEDEF)

Conservative temperature tracer.
"""
struct Temperature <: AbstractModelVariable end

field(::Temperature, model) = model.tracers.T
short_name(::Temperature) = "T"
standard_name(::Temperature) = "sea_water_conservative_temperature"
long_name(::Temperature) = "Conservative temperature"
units(::Temperature) = "°C"
color_map(::Temperature) = :thermal
symmetrize_limits(::Temperature) = false
function additional_attributes(::Temperature)
    Dict("unit_metadata" => "temperature: difference")
end

"""
$(TYPEDEF)

Salinity tracer.
"""
struct Salinity <: AbstractModelVariable end

field(::Salinity, model) = model.tracers.S
short_name(::Salinity) = "S"
standard_name(::Salinity) = "sea_water_salinity"
long_name(::Salinity) = "Salinity"
units(::Salinity) = "10⁻³"
color_map(::Salinity) = :haline
symmetrize_limits(::Salinity) = false

"""
$(TYPEDEF)

Specific kinetic energy.

## Details

Kinetic energy per unit mass, defined as

```math
eₖ(x, y, z, t) = \\frac{1}{2} (u^2 + v^2 + w^2)(x, y, z, t)
```
"""
struct SpecificKineticEnergy <: AbstractSpecificKineticEnergyVariable end

field(::SpecificKineticEnergy, model) = Field(sum(v^2 / 2 for v in model.velocities))
short_name(::SpecificKineticEnergy) = "eₖ"
standard_name(::SpecificKineticEnergy) = "specific_kinetic_energy_of_sea_water"
long_name(::SpecificKineticEnergy) = "Specific kinetic energy"

"""
$(TYPEDEF)

Specific turbulent kinetic energy.

## Details

Turbulent kinetic energy per unit mass, used as a tracer in CATKE turbulence closure schemes.
"""
struct SpecificTurbulentKineticEnergy <: AbstractSpecificKineticEnergyVariable end

field(::SpecificTurbulentKineticEnergy, model) = model.tracers.e
short_name(::SpecificTurbulentKineticEnergy) = "e"
standard_name(::SpecificTurbulentKineticEnergy) = "specific_turbulent_kinetic_energy_of_sea_water"
long_name(::SpecificTurbulentKineticEnergy) = "Specific turbulent kinetic energy"

"""
$(TYPEDEF)

Meridional overturning circulation (MOC) stream function.
    
## Details

Latitude-depth or y-depth fields corresponding to stream function of meridional overturning
circulation - computed here as vertically accumulated - that is cumulative vertical
integral with respect to depth - of zonally integrated meridional velocity component:

```math
\\Psi^M(\\varphi, z, t) = 
\\int_{0}^z \\int_{\\lambda_W}^{\\lambda_E} 
  v(\\lambda, \\varphi, z', t)
\\,\\mathrm{d}\\lambda \\,\\mathrm{d}z'
```

or in rectilinear coordinates:

```math
\\Psi^M(y, z, t) = 
\\int_{0}^z \\int_{x_W}^{x_E} 
  v(x, y, z', t)
\\,\\mathrm{d}x \\,\\mathrm{d}z'
```

The outputted field is scaled to be in sverdrup (10⁶ m³ s⁻¹) units.
"""
struct MOCStreamFunction <: AbstractVelocityStreamFunctionVariable end

function field(::MOCStreamFunction, model)
    # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
    # units 10⁶ m³ s⁻¹ = sverdrup
    Field(
        CumulativeIntegral(
            Field(Integral(1e-6 * model.velocities.v; dims=1)); dims=3, reverse=true
        ),
    )
end
short_name(::MOCStreamFunction) = "Ψᴹ"
standard_name(::MOCStreamFunction) = "ocean_meridional_overturning_streamfunction"
long_name(::MOCStreamFunction) = "MOC stream function"
units(::MOCStreamFunction) = "sverdrup"
function spatial_dimensions(::MOCStreamFunction)
    SpatialDimensions(
        ReducedSpatialDimension(), FreeSpatialDimension(), FreeSpatialDimension()
    )
end

"""
$(TYPEDEF)

Barotropic stream function.
    
## Details

Longitude-latitude or x-y fields corresponding to stream function of barotropic
velocity - computed here as zonally accumulated - that is cumulative integral with
respect to longitude - of depth integrated meridional velocity component:

```math
\\Psi^B(\\lambda, \\varphi, t) = 
\\int_{\\lambda_W}^{\\lambda}  \\int_{z_B}^{z_S} 
  v(\\lambda', \\varphi, z, t)
\\,\\mathrm{d}z \\,\\mathrm{d}\\lambda'
```

or in rectilinear coordinates:

```math
\\Psi^B(x, y, t) = 
\\int_{x_W}^{x}  \\int_{z_B}^{z_S} 
  v(x', y, z, t)
\\,\\mathrm{d}z \\,\\mathrm{d}x'
```

The outputted field is scaled to be in sverdrup (10⁶ m³ s⁻¹) units.
"""
struct BarotropicStreamFunction <: AbstractVelocityStreamFunctionVariable end

function field(::BarotropicStreamFunction, model)
    # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
    # units 10⁶ m³ s⁻¹ = sverdrup
    Field(CumulativeIntegral(Field(Integral(1e-6 * model.velocities.v; dims=3)); dims=1))
end
short_name(::BarotropicStreamFunction) = "Ψᴮ"
standard_name(::BarotropicStreamFunction) = "ocean_barotropic_streamfunction"
long_name(::BarotropicStreamFunction) = "Barotropic stream function"
units(::BarotropicStreamFunction) = "sverdrup"
function spatial_dimensions(::BarotropicStreamFunction)
    SpatialDimensions(
        FreeSpatialDimension(), FreeSpatialDimension(), ReducedSpatialDimension()
    )
end

"""
$(TYPEDEF)

Northward (meridional) heat transport.
    
$(TYPEDSIGNATURES)
    
## Details

Scalar field corresponding to northern (meridional) heat transport at a
given y / latitude coordinate.

The meridional heat transport is computed here as:

```math
Q(\\varphi, t) = 
\\int_{z_B}^{z_T} \\int_{\\lambda_W}^{\\lambda_E} 
  c_p \\rho v(\\lambda, \\varphi, z, t) T(\\lambda, \\varphi, z, t)
\\,\\mathrm{d}\\lambda \\,\\mathrm{d}z
```

or in rectilinear coordinates:

```math
Q(y, t) = 
\\int_{z_B}^{z_T} \\int_{x_W}^{x_E} 
  c_p \\rho v(x, y, z, t) T(x, y, z, t)
\\,\\mathrm{d}x \\,\\mathrm{d}z
```

where ``c_p`` and ``\\rho`` are reference values of the specific heat capacity
and density of sea water respectively.

The outputted heat transport is scaled to be in terrawatt (10¹² W) units.

$(TYPEDFIELDS)
"""
@kwdef struct NorthwardHeatTransport{T} <: AbstractModelVariable
    "Sea water reference density / kg m⁻³"
    sea_water_density::T = 1026.0
    "Sea water specific heat capacity / J K⁻¹ kg⁻¹"
    sea_water_specific_heat_capacity::T = 3991.0
end

function field(variable::NorthwardHeatTransport, model)
    Field(
        Integral(
            1e-12 * # Scale so integrated field is in units TW
            variable.sea_water_density *
            variable.sea_water_specific_heat_capacity *
            model.velocities.v *
            model.tracers.T;
            dims=(1, 3),
        ),
    )
end
short_name(::NorthwardHeatTransport) = "Q"
standard_name(::NorthwardHeatTransport) = "northward_ocean_heat_transport"
long_name(::NorthwardHeatTransport) = "Northward heat transport"
units(::NorthwardHeatTransport) = "TW"
color_map(::NorthwardHeatTransport) = :balance
function spatial_dimensions(::NorthwardHeatTransport)
    SpatialDimensions(
        ReducedSpatialDimension(), FreeSpatialDimension(), ReducedSpatialDimension()
    )
end

"""
$(TYPEDEF)

Eastward (zonal or x) component of barotropic velocity.
"""
struct EastwardBarotropicVelocity <: AbstractBarotropicVelocityVariable end

function field(::EastwardBarotropicVelocity, model)
    model.free_surface.barotropic_velocities.U
end
short_name(::EastwardBarotropicVelocity) = "U"
function standard_name(::EastwardBarotropicVelocity)
    "barotropic_eastward_sea_water_velocity"
end
long_name(::EastwardBarotropicVelocity) = "Barotropic eastward velocity"

"""
$(TYPEDEF)

Northward (meridional or y) component of barotropic velocity.
"""
struct NorthwardBarotropicVelocity <: AbstractBarotropicVelocityVariable end

function field(::NorthwardBarotropicVelocity, model)
    model.free_surface.barotropic_velocities.V
end
short_name(::NorthwardBarotropicVelocity) = "V"
function standard_name(::NorthwardBarotropicVelocity)
    "barotropic_northward_sea_water_velocity"
end
long_name(::NorthwardBarotropicVelocity) = "Barotropic northward velocity"

"""
$(TYPEDEF)

Free surface displacement.
"""
struct FreeSurfaceDisplacement <: AbstractModelVariable end

field(::FreeSurfaceDisplacement, model) = model.free_surface.displacement
short_name(::FreeSurfaceDisplacement) = "displacement"
standard_name(::FreeSurfaceDisplacement) = "sea_surface_height"
long_name(::FreeSurfaceDisplacement) = "Free surface displacement"
units(::FreeSurfaceDisplacement) = "m"
color_map(::FreeSurfaceDisplacement) = :balance
function spatial_dimensions(::FreeSurfaceDisplacement)
    SpatialDimensions(
        FreeSpatialDimension(), FreeSpatialDimension(), ReducedSpatialDimension()
    )
end

const VELOCITY_VARIABLES = (
    EastwardVelocity(), NorthwardVelocity(), UpwardVelocity()
)
const TRACER_VARIABLES = (Salinity(), Temperature())
const VELOCITY_AND_TRACER_VARIABLES = (VELOCITY_VARIABLES..., TRACER_VARIABLES...)
const BAROTROPIC_VELOCITY_VARIABLES = (
    EastwardBarotropicVelocity(), NorthwardBarotropicVelocity()
)

"""
$(SIGNATURES)

Per-variable output attributes for `variables` to save in NetCDF output.
"""
function output_attributes(variables::Vector{<:AbstractModelVariable})
    Dict(
        short_name(variable) => merge(
            Dict(
                "standard_name" => standard_name(variable),
                "long_name" => long_name(variable),
                "units" => units(variable),
            ),
            additional_attributes(variable),
        ) for variable in variables
    )
end

"""
$(SIGNATURES)

Automatic variable value limits for field for variable `variable` with 
data in `field_time_series`.

## Details

Computes extrema of field data and sets limits to these extrema, 
symmetrizing the upper and lower limits if appropriate for the variable.
"""
function auto_variable_limits(variable::AbstractModelVariable, field_timeseries)
    limits = extrema(interior(field_timeseries))
    symmetrize_limits(variable) ? (-maximum(abs.(limits)), maximum(abs.(limits))) : limits
end