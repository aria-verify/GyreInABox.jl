"""
$(TYPEDEF)

Parameters for baroclinic wind and buoyancy forced double-gyre ocean circulation model.
    
$(TYPEDSIGNATURES)

$(TYPEDFIELDS)
"""
@kwdef struct DoubleGyreParameters{T} <: AbstractParameters{T}
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
    "Average zonal wind velocity 10 meters above the ocean / m s⁻¹"
    u_10::T = 10.0
    "Dimensionless drag coefficient"
    c_d::T = 2e-3
    "Approximate average density of air at sea-level, kg m⁻³"
    ρ_a::T = 1.2
    "Latitude offset for zonal wind stress variation / °"
    φ_u::T = 15.0
    "Latitude scale for zonal wind stress variation / °"
    Lφ_u::T = 60.0
    "Latitude offset for reference salinity variation / °"
    φ_S::T = 15.0
    "Latitude scale for reference salinity variation / °"
    Lφ_S::T = 60.0
    "Bottom drag damping / m s⁻¹"
    μ::T = 1e-3
    "Sea water density / kg m⁻³"
    ρ_s::T = 1026.0
    "Sea water specific heat capacity / J K⁻¹ kg⁻¹"
    c_s::T = 3991.0
    "Reference polar surface temperature / °C"
    T_polar::T = 0.0
    "Reference equatorial surface temperature / °C"
    T_equatorial::T = 30.0
    "Reference abyssal ocean temperature / °C"
    T_abyssal::T = 2.0
    "Thermocline reference depth / m"
    z_thermocline::T = -1000.0
    "Thermocline reference length scale / m"
    ℓ_thermocline::T = 500.0
    "Halocline reference depth / m"
    z_halocline::T = -1000.0
    "Halocline reference length scale / m"
    ℓ_halocline::T = 500.0
    "Salinity reference level / g kg⁻¹"
    S_0::T = 35.0
    "Salinity latitude variation amplitude / g kg⁻¹"
    ΔS::T = 3.0
    "Salinity restoring timescale / s"
    τ_S::T = 90days
    "Temperature restoring timescale / s"
    τ_T::T = 30days
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
    T_surface = reference_surface_temperature(φ, p::DoubleGyreParameters)
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
function initial_salinity(φ, z, p::DoubleGyreParameters)
    S_surface = reference_surface_salinity(φ, p)
    sigmoidal_depth_profile(z, p.z_halocline, p.ℓ_halocline, p.S_0, S_surface)
end

# Use discrete form for field-dependent boundary condition due
# to https://github.com/CliMA/Oceananigans.jl/issues/4659

"""
Reference surface salinity for restoring surface salinity boundary condition as
function of latitude `φ` and parameters `p`.

$(SIGNATURES)
"""
@inline function reference_surface_salinity(φ, p::DoubleGyreParameters)
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
    i, j, grid, clock, model_fields, parameters_and_Δz::Tuple{DoubleGyreParameters{T},T}
) where {T}
    p, Δz = parameters_and_Δz
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds (model_fields.S[i, j, grid.Nz] - reference_surface_salinity(φ, p)) * Δz /
        p.τ_S
end

"""
Bottom surface drag on zonal velocity component in m² s⁻².

$(SIGNATURES)

## Details

Computes bottom surface drag at horizontal grid indices `i` and `j` for grid `grid` and
model clock `clock`, with current model fields `model_fields` and parameters `p`.
"""
@inline function bottom_zonal_drag(i, j, grid, clock, model_fields, p::DoubleGyreParameters)
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
    i, j, grid, clock, model_fields, p::DoubleGyreParameters
)
    @inbounds -p.μ * model_fields.v[i, j, 1]
end

"""
Reference surface temperature for restoring surface temperature boundary condition as
function of latitude `φ`` and parameters `p`.

$(SIGNATURES)
"""
@inline function reference_surface_temperature(φ, p::DoubleGyreParameters)
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
    i, j, grid, clock, model_fields, parameters_and_Δz::Tuple{DoubleGyreParameters{T},T}
) where {T}
    p, Δz = parameters_and_Δz
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    @inbounds (model_fields.T[i, j, grid.Nz] - reference_surface_temperature(φ, p)) * Δz /
        p.τ_T
end

"""
Maximal surface wind stress given parameters `p` in m² s⁻².

$(SIGNATURES)
"""
@inline function max_surface_wind_stress(p::DoubleGyreParameters)
    p.ρ_a / p.ρ_s * p.c_d * p.u_10^2
end

"""
Zonal wind stress applied as surface (flux) boundary condition to velocity field.

$(SIGNATURES)

## Details

Computes wind stress in zonal direction at longitude `λ` and latitude `φ` (both in
degrees), time `t` and with model parameters `p`.
"""
function zonal_wind_stress(λ, φ, t, p::DoubleGyreParameters)
    return max_surface_wind_stress(p) * cos(2π * (φ - p.φ_u) / p.Lφ_u)
end

function initialize!(model::Oceananigans.AbstractModel, parameters::DoubleGyreParameters)
    # Temperature initial condition: latitude and depth dependent thermocline profile
    T_initial(λ, φ, z) = initial_temperature(φ, z, parameters)
    # Salinity initial condition: latitude and depth dependent halocline profile
    S_initial(λ, φ, z) = initial_salinity(φ, z, parameters)
    set!(model; u=0.0, v=0.0, T=T_initial, S=S_initial)
    nothing
end

function grid(
    parameters::DoubleGyreParameters{T}, architecture::Oceananigans.AbstractArchitecture
) where {T}
    LatitudeLongitudeGrid(
        architecture;
        size=parameters.grid_size,
        longitude=parameters.longitude_interval,
        latitude=parameters.latitude_interval,
        z=k -> hyperbolically_spaced_faces(
            k,
            parameters.grid_size[3],
            parameters.depth_interval...,
            parameters.depth_stretching_factor,
        ),
        halo=parameters.halo_size,
        topology=(Bounded, Bounded, Bounded),
    )
end

function boundary_conditions(parameters::DoubleGyreParameters)
    Δz = minimum_zspacing(grid(parameters, CPU()))
    no_slip_bc = ValueBoundaryCondition(0)
    u_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(zonal_wind_stress; parameters=parameters),
        bottom=FluxBoundaryCondition(bottom_zonal_drag; discrete_form=true, parameters),
        north=no_slip_bc,
        south=no_slip_bc,
    )
    v_bcs = FieldBoundaryConditions(;
        bottom=FluxBoundaryCondition(
            bottom_meridional_drag; discrete_form=true, parameters
        ),
        east=no_slip_bc,
        west=no_slip_bc,
    )
    T_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(
            surface_temperature_flux; discrete_form=true, parameters=(parameters, Δz)
        ),
    )
    S_bcs = FieldBoundaryConditions(;
        top=FluxBoundaryCondition(
            surface_salinity_flux; discrete_form=true, parameters=(parameters, Δz)
        ),
    )
    return (; u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)
end

function buoyancy(::DoubleGyreParameters)
    SeawaterBuoyancy(; equation_of_state=TEOS10EquationOfState())
end

tracers(::DoubleGyreParameters) = (:T, :S, :e)

forcing(::DoubleGyreParameters) = NamedTuple()

momentum_advection(::DoubleGyreParameters) = Oceananigans.WENOVectorInvariant()

coriolis(::DoubleGyreParameters) = HydrostaticSphericalCoriolis()

closure(::DoubleGyreParameters) = CATKEVerticalDiffusivity()