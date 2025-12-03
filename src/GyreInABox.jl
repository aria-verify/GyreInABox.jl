"""
Oceananigans based model of an ocean gyre in a bounded domain.

## Details

Model of wind and buoyancy forced ocean gyre adapted from MITgcm
[baroclinic ocean gyre example from documentation](https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html)
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
using Oceananigans.Advection
using Oceananigans.Coriolis
using Oceananigans.Grids
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState
using CairoMakie
using Printf

export GyreInABoxParameters, GyreInABoxConfiguration
export HorizontalSlice, LongitudeDepthSlice, LatitudeDepthSlice, DepthTimeAveraged
export setup_model, initialize_model!, setup_simulation
export run_simulation, record_animation


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


abstract type AbstractOutputType{T} end

abstract type AbstractHorizontalOutputType{T} <: AbstractOutputType{T} end

abstract type AbstractVerticalOutputType{T} <: AbstractOutputType{T} end

"""
$(TYPEDEF)

Horizontal (latitude-longitude) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal slices through model velocity and tracer fields at specified depth
and interval.

$(TYPEDFIELDS)
"""
@kwdef struct HorizontalSlice{T} <: AbstractHorizontalOutputType{T}
    "Depth of slice / m"
    depth::T = 0.
    "Time interval to record output at / s"
    interval::T = 1day
end

"""
$(TYPEDEF)

Vertical (longitude-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified latitude
and interval.

$(TYPEDFIELDS)
"""
@kwdef struct LongitudeDepthSlice{T} <: AbstractVerticalOutputType{T}
    "Latitude of slice / °"
    latitude::T = 45.
    "Time interval to record output at / s"
    interval::T = 1day
end

"""
$(TYPEDEF)

Vertical (latitude-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified longitude
and interval.

$(TYPEDFIELDS)
"""
@kwdef struct LatitudeDepthSlice{T} <: AbstractVerticalOutputType{T}
    "Longitude of slice / °"
    longitude::T = 30.
    "Time interval to record output at / s"
    interval::T = 1day
end

"""
$(TYPEDEF)

Depth and time averaged output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal fields corresponding to depth averaged model velocity and traced 
fields temporally averaged over a specified window and a given temporal interval.

$(TYPEDFIELDS)
"""
@kwdef struct DepthTimeAveraged{T} <: AbstractHorizontalOutputType{T}
    "Time interval to record output at / s"
    interval::T = 30day
    "Time window to accumulate output averages over / s"
    average_window::T = 30day
end

SliceOutputType = Union{HorizontalSlice, LongitudeDepthSlice, LatitudeDepthSlice}

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
    MA<:Oceananigans.Advection.AbstractAdvectionScheme,
    C<:Oceananigans.Coriolis.AbstractRotation,
    TC<:Oceananigans.TurbulenceClosures.AbstractTurbulenceClosure
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
        DepthTimeAveraged()
    )
end

# Use discrete form for field-dependent boundary condition due
# to https://github.com/CliMA/Oceananigans.jl/issues/4659

"""
Reference surface salinity for restoring surface salinity boundary condition as
function of latitude `φ`` and parameters `p`.

$(SIGNATURES)
"""
@inline function reference_surface_salinity(φ, p::GyreInABoxParameters)
    p.S_0 + p.ΔS * sin(2π * (φ - p.φ_S) / p.Lφ_S)
end

"""
Surface salinity flux in g kg⁻¹ m s⁻¹.

$(SIGNATURES)

## Details

Relaxation boundary condition that restores surface salinity to latitude dependent
profile determined by `reference_surface_salinity`.
"""
@inline function surface_salinity_flux(
    i, j, grid, clock, model_fields, parameters_and_Δz
) 
    p, Δz = parameters_and_Δz
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds  (
        model_fields.S[i, j, grid.Nz] - reference_surface_salinity(φ, p)
    ) * Δz / p.τ_S
end

"""
Bottom surface drag on zonal velocity component in m² s⁻².

$(SIGNATURES)

## Details

Computes bottom surface drag at horizontal grid indices `i` and `j` for grid `grid` and
model clock `clock`, with current model fields `model_fields` and parameters `p`.
"""
@inline function bottom_zonal_drag(
    i, j, grid, clock, model_fields, p::GyreInABoxParameters
) 
    @inbounds -p.μ * model_fields.u[i, j, 1]
end

"""
Bottom surface drag on meridional velocity component in m² s⁻².

$(SIGNATURES)

## Details

Computes bottom surface drag at horizontal grid indices `i` and `j` for grid `grid` and
model clock `clock`, with current model fields `model_fields` and parameters `p`.
"""
@inline function bottom_meridional_drag(
    i, j, grid, clock, model_fields, p::GyreInABoxParameters
)
    @inbounds -p.μ * model_fields.v[i, j, 1]
end

"""
Reference surface temperature for restoring surface temperature boundary condition as
function of latitude `φ`` and parameters `p`.

$(SIGNATURES)
"""
@inline function reference_surface_temperature(φ, p::GyreInABoxParameters)
    p.T_polar + (p.T_equatorial - p.T_polar) * cosd(φ)^2
end

"""
Surface temperature flux in  K m s⁻¹.

$(SIGNATURES)

## Details

Relaxation boundary condition that restores surface temperature to latitude dependent
profile determined by `reference_surface_temperature`.
"""
@inline function surface_temperature_flux(
    i, j, grid, clock, model_fields, parameters_and_Δz
) 
    p, Δz = parameters_and_Δz
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds  (
        model_fields.T[i, j, grid.Nz] - reference_surface_temperature(φ, p)
    ) * Δz / p.τ_T
end

"""
Maximal surface wind stress given parameters `p` in m² s⁻².

$(SIGNATURES)
"""
@inline function max_surface_wind_stress(p::GyreInABoxParameters)
    p.ρ_a / p.ρ_s * p.c_d * p.u_10^2
end

"""
Zonal wind stress applied as surface (flux) boundary condition to velocity field.

$(SIGNATURES)

## Details

Computes wind stress in zonal direction at longitude `λ` and latitude `φ` (both in
degrees), time `t` and with model parameters `p`.
"""
function zonal_wind_stress(λ, φ, t, p::GyreInABoxParameters)
    return max_surface_wind_stress(p) * cos(2π * (φ - p.φ_u) / p.Lφ_u)
end

"""
Construct named tuple of boundary conditions for model given parameters `parameters`.

$(SIGNATURES)
"""
function boundary_conditions(
    parameters::GyreInABoxParameters{T}, grid::AbstractGrid
) where {T}
    Δz = minimum_zspacing(grid)
    no_slip_bc = ValueBoundaryCondition(0)
    u_bcs = FieldBoundaryConditions(
        top=FluxBoundaryCondition(zonal_wind_stress; parameters=parameters),
        bottom=FluxBoundaryCondition(bottom_zonal_drag; discrete_form=true, parameters),
        north=no_slip_bc,
        south=no_slip_bc
    )
    v_bcs = FieldBoundaryConditions(
        bottom=FluxBoundaryCondition(
            bottom_meridional_drag; discrete_form=true, parameters
        ),
        east=no_slip_bc,
        west=no_slip_bc
    )
    T_bcs = FieldBoundaryConditions(
        top=FluxBoundaryCondition(
            surface_temperature_flux; discrete_form=true, parameters=(parameters, Δz)
        )
    )
    S_bcs = FieldBoundaryConditions(
        top=FluxBoundaryCondition(
            surface_salinity_flux, discrete_form=true, parameters=(parameters, Δz)
        )
    )
    return (; u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)
end

"""
Compute hyperbolically spaced depth grid face coordinate at index `k` given model
configuration `configuration`.

$(SIGNATURES)
"""
function hyperbolically_spaced_faces(k::Int, configuration::GyreInABoxConfiguration)
    configuration.depth_interval[2] + (
        configuration.depth_interval[1] - configuration.depth_interval[2]
    ) * (
        1 - tanh(
            configuration.depth_stretching_factor * (k - 1) / configuration.grid_size[3]
        ) / tanh(configuration.depth_stretching_factor)
    )
end

"""
Set up ocean gyre model grid with configuration `configuration` optionally overriding
architecture in `configuration` with `architecture`.

$(SIGNATURES)
"""
function setup_grid(
    configuration::GyreInABoxConfiguration{T}; architecture = nothing
) where {T}
    LatitudeLongitudeGrid(
        isnothing(architecture) ? configuration.architecture : architecture,
        size=configuration.grid_size,
        longitude=configuration.longitude_interval,
        latitude=configuration.latitude_interval,
        z=k -> hyperbolically_spaced_faces(k, configuration),
        halo=configuration.halo_size,
        topology=(Bounded, Bounded, Bounded)
    )
end

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
Sigmoidal depth profile which smoothly varies from `f_top` at `z = 0` to `f_bottom` at
`z = -depth` with length scale for smooth transition `ℓ` and offset for midpoint of
smooth transition `z_0`.

$(SIGNATURES)
"""
function sigmoidal_depth_profile(z, z_0, ℓ, f_bottom, f_top)
    f_bottom + (f_top - f_bottom) * (1  + tanh((z - z_0) / ℓ)) / 2
end

"""
Initial temperature at latitude `φ` and depth `z` with parameters `p`.

$(SIGNATURES)

## Details

Computes a latitude and depth dependent initial temperature with thermocline like
profile which smoothly varies from a latitude dependent surface temperature to a
constant deep ocean temperature with a sigmoidal profile.
"""
function initial_temperature(φ, z, p)
    T_surface = reference_surface_temperature(φ, p)
    sigmoidal_depth_profile(z, p.z_thermocline, p.ℓ_thermocline, p.T_abyssal, T_surface)
end

"""
Initial salinity at latitude `φ` and depth `z` with parameters `p`.

$(SIGNATURES)

## Details

Computes a latitude and depth dependent initial salinity with halocline like profile
which smoothly varies from a latitude dependent surface salinity to a constant
deep ocean salinity with a sigmoidal profile.
"""
function initial_salinity(φ, z, p)
    S_surface = reference_surface_salinity(φ, p)
    sigmoidal_depth_profile(z, p.z_halocline, p.ℓ_halocline, p.S_0, S_surface)
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


"""
Symbol label for output type to use in naming output file and registering output writer.

$(TYPEDSIGNATURES)
"""
function label(::AbstractOutputType) end

label(::HorizontalSlice) = :horizontal_slice
label(::LongitudeDepthSlice) = :longitude_depth_slice
label(::LatitudeDepthSlice) = :latitude_depth_slice
label(::DepthTimeAveraged) = :depth_time_averaged

"""
Spatial grid indices output type records fields at.

$(TYPEDSIGNATURES)
"""
function indices(type::AbstractOutputType, grid) end

indices(type::HorizontalSlice, grid) = (
    :, :, clamp(searchsortedfirst(znodes(grid, Face()), type.depth), 1:grid.Nz)
)
indices(type::LongitudeDepthSlice, grid) = (
    :, clamp(searchsortedfirst(φnodes(grid, Face()), type.latitude), 1:grid.Ny), :
)
indices(type::LatitudeDepthSlice, grid) = (
    clamp(searchsortedfirst(λnodes(grid, Face()), type.longitude), 1:grid.Nx), :, :
)
indices(::DepthTimeAveraged, grid) = (:, :, :)

"""
Named tuple of output variables (fields) to record for output type.
    
$(TYPEDSIGNATURES)
"""
function outputs(::AbstractOutputType, model) end

outputs(::SliceOutputType, model) = merge(model.velocities, model.tracers)
outputs(::DepthTimeAveraged, model) = NamedTuple(
    name => Field(Average(variable, dims=3))
    for (name, variable) in pairs(merge(model.velocities, model.tracers))
)

"""
Time schedule to record output type at.

$(TYPEDSIGNATURES)
"""
function schedule(::AbstractOutputType) end

schedule(type::SliceOutputType) = TimeInterval(type.interval)
schedule(type::DepthTimeAveraged) = AveragedTimeInterval(
    type.interval, window=type.average_window
)

"""
Filename to record outputs to.

$(TYPEDSIGNATURES)

## Details

Output filename stem is specified in `configuration`. For an output type `type`
a label computed using `label` function is appended on to stem and file extension is
specified by `extension` added.
"""
output_filename(
    configuration::GyreInABoxConfiguration, label::Symbol, extension::String
) = "$(configuration.output_filename)_$(label).$(extension)"
output_filename(
    configuration::GyreInABoxConfiguration,
    type::AbstractOutputType,
    extension::String="jld2"
) = output_filename(configuration, label(type), extension)

"""
Register output writers for output types specified in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_output_writers!(
    simulation::Oceananigans.Simulation,
    configuration::GyreInABoxConfiguration{T}
) where {T}

    model = simulation.model
    # Create a grid on CPU if not already to avoid issues with computing output field
    # indices from grid using scalar operations on GPU grids
    grid = (
        isa(configuration.architecture, CPU) 
        ? model.grid 
        : setup_grid(configuration, architecture=CPU())
    )

    for output_type in configuration.output_types
        if output_type.interval <= configuration.simulation_time
            simulation.output_writers[label(output_type)] = JLD2Writer(
                model,
                outputs(output_type, model),
                filename=output_filename(configuration, output_type),
                indices=indices(output_type, grid),
                schedule=schedule(output_type),
                overwrite_existing=true,
                with_halos=true,
            )
        end
    end

    nothing
end

"""
Add callback for progress updates as configured in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_progress_message_callback!(
    simulation::Oceananigans.Simulation,
    configuration::GyreInABoxConfiguration{T}
) where {T}
    message_string_format = Printf.Format(
        "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n  " *
        "max(|u|) = %.1e ms⁻¹, max(|v|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, " *
        "max(|S|) = %.3e g/kg, max(|T|) = %.3e ᵒC, max(|e|) = %.3e m²s⁻²"
    )
    progress_message(sim) = println(
        Printf.format(
            message_string_format,
            iteration(sim),
            prettytime(sim),
            prettytime(sim.Δt),
            prettytime(sim.run_wall_time),
            maximum(abs, sim.model.velocities.u),
            maximum(abs, sim.model.velocities.v),
            maximum(abs, sim.model.velocities.w),
            maximum(abs, sim.model.tracers.S),
            maximum(abs, sim.model.tracers.T),
            maximum(abs, sim.model.tracers.e)
        )
    )
    add_callback!(
        simulation,
        progress_message,
        IterationInterval(configuration.progress_message_interval)
    )
    nothing
end

"""
Set up simulation for model `model` with configuration `configuration`.
    
$(SIGNATURES)
"""
function setup_simulation(
    model::Oceananigans.AbstractModel,
    configuration::GyreInABoxConfiguration{T}
) where {T}
    simulation = Simulation(
        model,
        Δt=configuration.initial_timestep,
        stop_time=configuration.simulation_time
    )
    wizard = TimeStepWizard(
        cfl=configuration.target_cfl,
        max_change=configuration.wizard_max_change,
        max_Δt=configuration.maximum_timestep
    )
    simulation.callbacks[:wizard] = Callback(
        wizard, IterationInterval(configuration.wizard_update_interval)
    )
    add_progress_message_callback!(simulation, configuration)
    add_output_writers!(simulation, configuration)
    simulation
end

"""
Setup and initialize model then setup and run simulation with parameters `parameters`
and configuration `configuration`.

$(SIGNATURES)
"""
function run_simulation(
    parameters::GyreInABoxParameters{T},
    configuration::GyreInABoxConfiguration{T}
) where {T}
    model = setup_model(parameters, configuration)
    initialize_model!(model, parameters)
    simulation = setup_simulation(model, configuration)
    run!(simulation)
end

"""
Label for horizontal axis for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_xlabel(::AbstractOutputType) end

axis_xlabel(::AbstractHorizontalOutputType) = "Longitude λ (ᵒ)"
axis_xlabel(::LongitudeDepthSlice) = "Longitude λ (ᵒ)"
axis_xlabel(::LatitudeDepthSlice) = "Latitude ϕ (ᵒ)"

"""
Label for vertical axis for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_ylabel(::AbstractOutputType) end

axis_ylabel(::AbstractHorizontalOutputType) = "Latitude ϕ (ᵒ)"
axis_ylabel(::AbstractVerticalOutputType) = "Depth z (m)"


"""
Aspect ratio for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
axis_aspect_ratio(configuration, output_type::AbstractOutputType) = AxisAspect(1)
axis_aspect_ratio(configuration, output_type::AbstractHorizontalOutputType) = AxisAspect(
    abs(configuration.longitude_interval[1] - configuration.longitude_interval[2]) /
    abs(configuration.latitude_interval[1] - configuration.latitude_interval[2])
)

"""
Axis limits for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_limits(configuration::GyreInABoxConfiguration, ::AbstractOutputType) end

axis_limits(configuration::GyreInABoxConfiguration, ::AbstractHorizontalOutputType) = (
    configuration.longitude_interval, configuration.latitude_interval
)
axis_limits(configuration::GyreInABoxConfiguration, ::LongitudeDepthSlice) = (
    configuration.longitude_interval, configuration.depth_interval
)
axis_limits(configuration::GyreInABoxConfiguration, ::LatitudeDepthSlice) = (
    configuration.latitude_interval, configuration.depth_interval
)

"""
Record an animation of fields recorded as model output.

$(SIGNATURES)

## Details

Generates and record to file an animation of model outputs for a model configuration
`configuration` and output type `output_type`. The simulation must have already been
run for this configuration and with specified output type active.

Animation is generated using CairoMakie with a figure window size in pixels of
`(figure_size)`, frame rate in frames per second of `frame_rate` and limits for the
zonal velocity, meridional velocity, vertical velocity, temperature, salinity and
turbulent kinetic energy field value color mappings specified by `u_limits`, `v_limits`,
`w_limits`, `T_limits`, `S_limits` and `e_limits` respectively.
"""
function record_animation(
    configuration::GyreInABoxConfiguration{T},
    output_type::AbstractOutputType;
    figure_size::Tuple = (1800, 900),
    frame_rate::Int = 8,
    u_limits::Tuple = (-2., 2.),
    v_limits::Tuple = (-2., 2.),
    w_limits::Tuple = (-1e-2, 1e-2),
    T_limits::Tuple = (2., 30.),
    S_limits::Tuple = (32., 38.),
    e_limits::Tuple = (0., 1e-2)
) where {T}
    filepath = output_filename(configuration, output_type)

    time_series = (
        u=FieldTimeSeries(filepath, "u"),
        v=FieldTimeSeries(filepath, "v"),
        w=FieldTimeSeries(filepath, "w"),
        T=FieldTimeSeries(filepath, "T"),
        S=FieldTimeSeries(filepath, "S"),
        e=FieldTimeSeries(filepath, "e")
    )

    times = time_series.u.times

    n = Observable(1)

    uₙ = @lift time_series.u[$n]
    vₙ = @lift time_series.v[$n]
    wₙ = @lift time_series.w[$n]
    Tₙ = @lift time_series.T[$n]
    Sₙ = @lift time_series.S[$n]
    eₙ = @lift time_series.e[$n]

    fig = Figure(size=figure_size)

    axis_kwargs = (
        xlabel=axis_xlabel(output_type),
        ylabel=axis_ylabel(output_type),
        aspect=axis_aspect_ratio(configuration, output_type),
        limits=axis_limits(configuration, output_type)
    )

    ax_u = Axis(fig[2, 1]; title="Zonal velocity u", axis_kwargs...)
    ax_v = Axis(fig[2, 3]; title="Meridional velocity v", axis_kwargs...)
    ax_w = Axis(fig[2, 5]; title="Vertical velocity w", axis_kwargs...)
    ax_T = Axis(fig[3, 1]; title="Temperature T", axis_kwargs...)
    ax_S = Axis(fig[3, 3]; title="Salinity S", axis_kwargs...)
    ax_e = Axis(fig[3, 5]; title="Turbulent kinetic energy e", axis_kwargs...)

    title = @lift @sprintf("t = %s", prettytime(times[$n]))

    hm_u = heatmap!(ax_u, uₙ; colormap=:balance, colorrange=u_limits)
    Colorbar(fig[2, 2], hm_u; label="m s⁻¹")

    hm_v = heatmap!(ax_v, vₙ; colormap=:balance, colorrange=v_limits)
    Colorbar(fig[2, 4], hm_v; label="m s⁻¹")
    
    hm_w = heatmap!(ax_w, wₙ; colormap=:balance, colorrange=w_limits)
    Colorbar(fig[2, 6], hm_w; label="m s⁻¹")

    hm_T = heatmap!(ax_T, Tₙ; colormap=:thermal, colorrange=T_limits)
    Colorbar(fig[3, 2], hm_T; label="ᵒC")

    hm_S = heatmap!(ax_S, Sₙ; colormap=:haline, colorrange=S_limits)
    Colorbar(fig[3, 4], hm_S; label="g / kg")
    
    hm_e = heatmap!(ax_e, eₙ; colormap=:thermal, colorrange=e_limits)
    Colorbar(fig[3, 6], hm_e; label="m² s⁻²")

    fig[1, 1:6] = Label(fig, title, fontsize=24, tellwidth=false)

    frames = 1:length(times)

    output_file = output_filename(configuration, output_type, "mp4")
    @info "Recording an animation of $(label(output_type)) to $(output_file)..."

    CairoMakie.record(fig, output_file, frames, framerate=frame_rate) do i
        (i % 10 == 0) && @printf "Plotting frame %i of %i\n" i frames[end]
        n[] = i
    end

    fig
end

end
