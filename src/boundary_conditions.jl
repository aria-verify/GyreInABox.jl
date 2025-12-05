# Use discrete form for field-dependent boundary condition due
# to https://github.com/CliMA/Oceananigans.jl/issues/4659

"""
Reference surface salinity for restoring surface salinity boundary condition as
function of latitude `φ` and parameters `p`.

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
