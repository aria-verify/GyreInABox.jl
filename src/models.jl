"""Parameters for an ocean gyre model.

Models should define a concrete subtype of this type.
"""
abstract type AbstractParameters{T} end

"""
    $(FUNCTIONNAME)(model::Oceananigans.AbstractModel, parameters)

Initialize state of ocean gyre model `model` with parameters `parameters`.
"""
function initialize! end

"""
    $(FUNCTIONNAME)(parameters, architecture::Oceananigans.AbstractArchitecture)

Set up ocean gyre model grid with parameters `parameters` on architecture `architecture`.
"""
function grid end

"""
    $(FUNCTIONNAME)(parameters)

Construct named tuple of boundary conditions for ocean gyre model given parameters `parameters`.
"""
function boundary_conditions end

"""
    $(FUNCTIONNAME)(parameters)

Construct named tuple of forcings for ocean gyre model given parameters `parameters`.
"""
function forcing end

"""
    $(FUNCTIONNAME)(parameters)

Construct buoyancy formulation for ocean gyre model given parameters `parameters`.
"""
function buoyancy end

"""
    $(FUNCTIONNAME)(parameters)

Construct turbulence closure(s) for ocean gyre model given parameters `parameters`.
"""
function closure end

"""
    $(FUNCTIONNAME)(parameters)

Construct Coriolis force formulation for ocean gyre model given parameters `parameters`.
"""
function coriolis end

"""
    $(FUNCTIONNAME)(parameters)

Construct tuple of tracer field names for ocean gyre model given parameters `parameters`.
"""
function tracers end

"""
    $(FUNCTIONNAME)(parameters)

Construct momentum advection scheme for ocean gyre model with parameters `parameters`.
"""
function momentum_advection end

"""
Set up ocean gyre model with parameters `parameters` on grid with architecture `architecture`.

$(SIGNATURES)
"""
function setup_model(parameters::AbstractParameters; architecture=CPU())
    HydrostaticFreeSurfaceModel(
        grid(parameters, architecture);
        momentum_advection=momentum_advection(parameters),
        coriolis=coriolis(parameters),
        closure=closure(parameters),
        buoyancy=buoyancy(parameters),
        boundary_conditions=boundary_conditions(parameters),
        forcing=forcing(parameters),
        tracers=tracers(parameters),
    )
end