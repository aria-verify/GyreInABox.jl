"""
$(TYPEDEF)

Parameters for simplified ocean gyre model from Spall (2011).
    
$(TYPEDSIGNATURES)
    
## Details

Real-valued parameters of model controlling initial and boundary conditions.

## References

1. Spall, M.A., 2011. On the Role of Eddies and Surface Forcing in the Heat Transport 
   and Overturning Circulation in Marginal Seas. Journal of Climate, 24, 4844--4858,
   https://doi.org/10.1175/2011JCLI4130.

$(TYPEDFIELDS)
"""
@kwdef struct Spall2011Parameters{T} <: AbstractParameters{T}
    "Grid dimensions in x, y and depth"
    grid_size::NTuple{3,Int} = (200, 400, 20)
    "Dimensions of grid halo region in x, y and depth"
    halo_size::NTuple{3,Int} = (7, 7, 4)
    "β-plane Coriolis offset parameters / s⁻¹"
    coriolis_offset::T = 1.2e-4
    "β-plane Coriolis coefficient parameter / m⁻¹s⁻¹"
    coriolis_coefficient::T = 2e-11
    "Zonal wind stress amplitue / N m⁻²"
    zonal_wind_stress::T = 0.15
    "Meridional wind stress amplitue / N m⁻²"
    meridional_wind_stress::T = 0.0
    "Surface temperature restoring strength / W m⁻² K⁻¹"
    surface_temperature_restoring_strength::T = 20.0
    "Buoyancy vertical stratification coefficient / s⁻²"
    vertical_stratification::T = 2e-6
    "Vertical scalar viscosity turbulence closure coefficient / m² s⁻¹"
    vertical_viscosity_coefficient::T = 1e-5
    "Vertical scalar diffusivity turbulence closure coefficient / m² s⁻¹"
    vertical_diffusivity_coefficient::T = 1e-5
    "Thermal expansion coefficient / kg m⁻³ K⁻¹"
    thermal_expansion_coefficient::T = 0.2
    "Sea water reference density / kg m⁻³"
    sea_water_density::T = 1026.0
    "Sea water specific heat capacity / J K⁻¹ kg⁻¹"
    sea_water_heat_capacity::T = 3991.0
    "Northern boundary surface temperature / °C"
    northern_surface_temperature::T = 2.0
    "Southern boundary surface temperature / °C"
    southern_surface_temperature::T = 10.0
    "Southern region temperature relaxation time scale / s"
    southern_region_temperature_relaxation_time::T = 20day
    "Southern region extent / m"
    southern_region_extent::T = 200kilometers
    "Southern region boundary smoothing window / m"
    southern_boundary_window_width::T = 0.0
    "Scale factor for exponentially spaced depth grid"
    depth_grid_scale_factor::T = 825.0
    "Domain size in x dimension / m"
    domain_size_x::T = 1000kilometers
    "Domain size in y dimension / m"
    domain_size_y::T = 2000kilometers
    "Depth of bottom of domain (most negative z) / m"
    bottom_depth::T = -2kilometers
    "Location of center of sill on sea floor along y dimension / m"
    sill_center_y::T = 1200kilometers
    "Width of sill on sea floor / m"
    sill_width::T = 400kilometers
    "Height of sill on sea floor / m"
    sill_height::T = 1kilometers
    "Width of slope on side walls of domain / m"
    side_slope_width::T = 140kilometers
    "Depth at which slope on side walls starts / m"
    side_slope_top_depth::T = 50meters
    "Width of slope on top wall of domain / m"
    top_slope_width::T = 20kilometers
    "Whether to use CATKE rather than scalar vertical diffusivity turbulence closure"
    use_catke_closure::Bool = true
    "Whether to initialise with reference surface temperature or constant"
    initialize_with_reference_surface_temperature::Bool = true
    "Surface wind forcing ramp-up time scale / s"
    surface_wind_forcing_ramp_up_timescale::T = 10day
end

function zonal_surface_wind_stress(x, y, t, parameters::Spall2011Parameters)
    (parameters.zonal_wind_stress / parameters.sea_water_density) *
    smooth_step(t / parameters.surface_wind_forcing_ramp_up_timescale) *
    cospi(y / parameters.domain_size_y)
end

function meridional_surface_wind_stress(x, y, t, parameters::Spall2011Parameters)
    (parameters.meridional_wind_stress / parameters.sea_water_density) *
    smooth_step(t / parameters.surface_wind_forcing_ramp_up_timescale) *
    cospi(x / parameters.domain_size_x)
end

function sill_profile(y, center, width, height)
    height * cospi((y - center) / width)^2
end

function side_wall_profile(x, width, height)
    x < width ? height * (1 - (x / width)) : 0
end

function two_basin_bathymetry(x, y, parameters::Spall2011Parameters)
    depth = parameters.bottom_depth
    northern_basin_radius = parameters.domain_size_x / 2
    sill_northern_limit = parameters.sill_center_y + parameters.sill_width / 2
    sill_southern_limit = parameters.sill_center_y - parameters.sill_width / 2
    if y > parameters.domain_size_y - northern_basin_radius
        # in northern basin rounded region
        x_basin = x - northern_basin_radius
        y_basin = y - (parameters.domain_size_y - northern_basin_radius)
        r_basin = sqrt(x_basin^2 + y_basin^2)
        if r_basin > northern_basin_radius
            depth = 0
        else
            depth += side_wall_profile(
                northern_basin_radius - r_basin,
                parameters.top_slope_width +
                (x_basin / r_basin)^2 *
                (parameters.side_slope_width - parameters.top_slope_width),
                parameters.side_slope_top_depth - parameters.bottom_depth,
            )
        end
    elseif y > sill_southern_limit && y < sill_northern_limit
        # in sill region
        sill_height = sill_profile.(
            y, parameters.sill_center_y, parameters.sill_width, parameters.sill_height
        )
        x_boundary = min(x, parameters.domain_size_x - x)
        side_slope_height = side_wall_profile(
            x_boundary,
            parameters.side_slope_width,
            parameters.side_slope_top_depth - parameters.bottom_depth,
        )
        depth += max(sill_height, side_slope_height)
    else
        # in southern basin or non-rounded northern basin region
        x_boundary = min(x, parameters.domain_size_x - x)
        if x_boundary < parameters.side_slope_width
            depth += side_wall_profile(
                x_boundary,
                parameters.side_slope_width,
                parameters.side_slope_top_depth - parameters.bottom_depth,
            )
        end
    end
    return depth
end

@inline function reference_surface_temperature(x, y, p::Spall2011Parameters)
    p.southern_surface_temperature +
    (max(y - p.southern_region_extent, 0) / (p.domain_size_y - p.southern_region_extent)) *
    (p.northern_surface_temperature - p.southern_surface_temperature)
end

@inline function surface_temperature_flux(
    i, j, grid, clock, model_fields, p::Spall2011Parameters
)
    x = xnode(i, j, 1, grid, Center(), Center(), Center())
    y = ynode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds (model_fields.T[i, j, grid.Nz] - reference_surface_temperature(x, y, p)) * (
        p.surface_temperature_restoring_strength /
        (p.sea_water_density * p.sea_water_heat_capacity)
    )
end

@inline function southern_region_mask(x, y, z, p::Spall2011Parameters)
    if p.southern_boundary_window_width == 0
        y < p.southern_region_extent
    else
        # Use a smooth step function at boundary of mask
        # Distance from boundary normalized by window width and offset by 0.5 as
        # smooth_step interopolates from 0 to 1 over (0, 1) interval
        d = (p.southern_region_extent - y) / p.southern_boundary_window_width + 0.5
        smooth_step(d)
    end
end

@inline function vertically_stratified_temperature(
    z, surface_temperature, p::Spall2011Parameters
)
    surface_temperature +
    z * p.sea_water_density * p.vertical_stratification /
    (p.thermal_expansion_coefficient * Oceananigans.defaults.gravitational_acceleration)
end

@inline function southern_region_temperature_target(x, y, z, t, p::Spall2011Parameters)
    vertically_stratified_temperature(z, p.southern_surface_temperature, p)
end

@inline function southern_region_temperature_forcing(
    i, j, k, grid, clock, model_fields, parameters::Spall2011Parameters
)
    t = clock.time
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    return @inbounds(
        southern_region_mask(x, y, z, parameters) * (
            southern_region_temperature_target(x, y, z, t, parameters) -
            model_fields.T[i, j, k]
        ) / parameters.southern_region_temperature_relaxation_time
    )
end

function boundary_conditions(parameters::Spall2011Parameters{T}) where {T}
    u_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(zonal_surface_wind_stress; parameters=parameters)
    )
    v_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(meridional_surface_wind_stress; parameters=parameters)
    )
    T_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(surface_temperature_flux; discrete_form=true, parameters)
    )
    return (; u=u_bcs, v=v_bcs, T=T_bcs)
end

function forcing(parameters::Spall2011Parameters)
    (; T=Forcing(southern_region_temperature_forcing; discrete_form=true, parameters))
end

function grid(
    parameters::Spall2011Parameters, architecture::Oceananigans.AbstractArchitecture
)
    underlying_grid = RectilinearGrid(
        architecture;
        size=parameters.grid_size,
        x=(0, parameters.domain_size_x),
        y=(0, parameters.domain_size_y),
        z=ExponentialDiscretization(
            parameters.grid_size[3],
            parameters.bottom_depth,
            0;
            scale=parameters.depth_grid_scale_factor,
        ),
        halo=parameters.halo_size,
        topology=(Bounded, Bounded, Bounded),
    )
    ImmersedBoundaryGrid(
        underlying_grid, GridFittedBottom((x, y) -> two_basin_bathymetry(x, y, parameters))
    )
end

function buoyancy(parameters::Spall2011Parameters)
    equation_of_state = LinearEquationOfState(;
        thermal_expansion=(
            parameters.thermal_expansion_coefficient / parameters.sea_water_density
        ),
        haline_contraction=0.0,
    )
    SeawaterBuoyancy(; equation_of_state)
end

function closure(parameters::Spall2011Parameters)
    if parameters.use_catke_closure
        CATKEVerticalDiffusivity()
    else
        VerticalScalarDiffusivity(;
            ν=parameters.vertical_viscosity_coefficient,
            κ=parameters.vertical_diffusivity_coefficient,
        )
    end
end

function coriolis(parameters::Spall2011Parameters)
    BetaPlane(parameters.coriolis_offset, parameters.coriolis_coefficient)
end
tracers(parameters::Spall2011Parameters) = (:T, :S)
momentum_advection(parameters::Spall2011Parameters) = Oceananigans.WENOVectorInvariant()

function initialize!(model::Oceananigans.AbstractModel, parameters::Spall2011Parameters)
    function T_initial(x, y, z)
        surface_temperature = if parameters.initialize_with_reference_surface_temperature
            reference_surface_temperature(x, y, parameters)
        else
            parameters.southern_surface_temperature, parameters
        end
        vertically_stratified_temperature(z, surface_temperature, parameters)
    end
    set!(model; u=0.0, v=0.0, T=T_initial, S=0.0)
    nothing
end
