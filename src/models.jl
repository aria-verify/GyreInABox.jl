"""
Set up ocean gyre model with parameters `parameters` and configuration `configuration`.

$(SIGNATURES)
"""
function setup_model(
    parameters::GyreInABoxParameters{T}, configuration::GyreInABoxConfiguration{T}
) where {T}
    grid = setup_grid(configuration)
    HydrostaticFreeSurfaceModel(
        ; grid,
        momentum_advection=configuration.momentum_advection,
        coriolis=configuration.coriolis,
        closure=configuration.closure,
        buoyancy=SeawaterBuoyancy(; equation_of_state=configuration.equation_of_state),
        boundary_conditions=boundary_conditions(parameters, grid),
        tracers=configuration.tracers
    )
end

"""
Initialize state of ocean gyre model `model` with parameters `parameters`.

$(SIGNATURES)
"""
function initialize_model!(
    model::Oceananigans.AbstractModel,
    parameters::GyreInABoxParameters{T}
) where {T}
    # Temperature initial condition: latitude and depth dependent thermocline profile
    T_initial(λ, φ, z) = initial_temperature(φ, z, parameters)
    # Salinity initial condition: latitude and depth dependent halocline profile
    S_initial(λ, φ, z) = initial_salinity(φ, z, parameters)
    set!(model, u=0., v=0., T=T_initial, S=S_initial)
    nothing
end
