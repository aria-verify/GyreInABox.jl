"""
$(TYPEDEF)

Parameters for ocean gyre model.
    
$(TYPEDSIGNATURES)
    
## Details

Real-valued parameters of model controlling initial and boundary conditions.

$(TYPEDFIELDS)
"""
@kwdef struct GyreInABoxParameters{T}
    "Average zonal wind velocity 10 meters above the ocean / m s⁻¹"
    u_10::T = 10.
    "Dimensionless drag coefficient"  
    c_d::T = 2e-3
    "Approximate average density of air at sea-level, kg m⁻³"
    ρ_a::T = 1.2
    "Latitude offset for zonal wind stress variation / °"
    φ_u::T = 15.
    "Latitude scale for zonal wind stress variation / °"  
    Lφ_u::T = 60.
    "Latitude offset for reference salinity variation / °"
    φ_S::T = 15.
    "Latitude scale for reference salinity variation / °"  
    Lφ_S::T = 60.
    "Bottom drag damping / m s⁻¹"
    μ::T = 1e-3
    "Sea water density / kg m⁻³"
    ρ_s::T = 1026.
    "Sea water specific heat capacity / J K⁻¹ kg⁻¹"
    c_s::T = 3991.
    "Reference polar surface temperature / °C"
    T_polar::T = 0.
    "Reference equatorial surface temperature / °C"
    T_equatorial::T = 30.
    "Reference abyssal ocean temperature / °C"
    T_abyssal::T = 2.
    "Thermocline reference depth / m"
    z_thermocline::T = -1000.
    "Thermocline reference length scale / m"
    ℓ_thermocline::T = 500.
    "Halocline reference depth / m"
    z_halocline::T = -1000.
    "Halocline reference length scale / m"
    ℓ_halocline::T = 500.
    "Salinity reference level / g kg⁻¹"
    S_0::T = 35.
    "Salinity latitude variation amplitude / g kg⁻¹"
    ΔS::T = 3.
    "Salinity restoring timescale / s"
    τ_S::T = 90days
    "Temperature restoring timescale / s"
    τ_T::T = 30days
end

"""
$(TYPEDEF)

Configuration for ocean gyre model.
    
$(TYPEDSIGNATURES)
    
## Details

Variables defining overall configuration of model such as spatial and temporal domains
and corresponding discretizations, numerical schemes to use for different model
components and outputs to record.

Compared to the `GyreInABoxParameters` the fields here mostly have discrete values and
would not be typically calibrated to data.

$(TYPEDFIELDS)
"""
@kwdef struct GyreInABoxConfiguration{
    T,
    A<:Oceananigans.AbstractArchitecture,
    EOS,
    MA<:AbstractAdvectionScheme,
    C<:AbstractRotation,
    TC<:AbstractTurbulenceClosure
}
    "Computational architecture to run simulation on"
    architecture::A = Oceananigans.CPU()
    "Momentum advection scheme"
    momentum_advection::MA = Oceananigans.WENOVectorInvariant()
    "Buoyancy equation of state"
    equation_of_state::EOS = TEOS10EquationOfState()
    "Coriolis force"
    coriolis::C = HydrostaticSphericalCoriolis()
    "Turbulence closure"
    closure::TC = CATKEVerticalDiffusivity()
    "Tracer variables - needs to match turbulence closure and buoyancy formulation."
    tracers::Tuple = (:T, :S, :e)
    "Grid dimensions in longitude, latitude and depth"
    grid_size::Tuple = (60, 60, 15)
    "Extent of spatial domain in longitude / °"
    longitude_interval::Tuple = (0, 60)
    "Extent of spatial domain in latitude / °"
    latitude_interval::Tuple = (15, 75)
    "Extent of spatial domain in depth / m"
    depth_interval::Tuple = (-1.8kilometers, 0)
    "Stretching factor for hyperbolic spaced depth grid"
    depth_stretching_factor::T = 1.2
    "Dimensions of grid halo region in longitude, latitude and depth"
    halo_size::Tuple = (6, 6, 3)
    "Time to simulate for / s"
    simulation_time::T = 60day
    "Initial time step to use / s"
    initial_timestep::T = 10minute
    "Maximum time step to use in time step adaptation / s"
    maximum_timestep::T = 30minute
    "Stem of output file names"
    output_filename::String = "gyre_model"
    "Iteration interval between progress messages"
    progress_message_interval::Int = 40
    "Target (advective) CFL number for time stepping wizard"
    target_cfl::T = 0.2
    "Update (iteration) interval for time stepping wizard"
    wizard_update_interval::Int = 10
    "Maximum relative time step change in each wizard update"
    wizard_max_change::T = 1.5
    "Output types to record during simulation"
    output_types::Tuple = (
        HorizontalSlice(),
        LongitudeDepthSlice(),
        LatitudeDepthSlice(),
        FreeSurfaceFields(),
        MOCStreamFunction()
    )
end
