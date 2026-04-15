abstract type AbstractModelOutput{S} end

abstract type AbstractHorizontalModelOutput{S} <: AbstractModelOutput{S} end

abstract type AbstractVerticalModelOutput{S} <: AbstractModelOutput{S} end

abstract type AbstractScalarModelOutput{S} <: AbstractModelOutput{S} end

"""
$(TYPEDEF)

Horizontal (latitude-longitude) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal slices through model velocity and tracer fields at specified depth.

$(TYPEDFIELDS)
"""
@kwdef struct HorizontalSlice{S,T} <: AbstractHorizontalModelOutput{S}
    "Depth of slice / m"
    depth::T = 0.0
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Vertical (longitude-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified latitude.

$(TYPEDFIELDS)
"""
@kwdef struct LongitudeDepthSlice{S,T} <: AbstractVerticalModelOutput{S}
    "Latitude of slice / °"
    latitude::T = 45.0
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Vertical (latitude-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified longitude.

$(TYPEDFIELDS)
"""
@kwdef struct LatitudeDepthSlice{S,T} <: AbstractVerticalModelOutput{S}
    "Longitude of slice / °"
    longitude::T = 30.0
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Vertical (x-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified y coordinate.

$(TYPEDFIELDS)
"""
@kwdef struct XDepthSlice{S,T} <: AbstractVerticalModelOutput{S}
    "y coordinate of slice / m"
    y::T = 1000.0
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Vertical (y-depth) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records vertical slices through model velocity and tracer fields at specified x coordinate.

$(TYPEDFIELDS)
"""
@kwdef struct YDepthSlice{S,T} <: AbstractVerticalModelOutput{S}
    "x coordinate of slice / m"
    x::T = 500.0
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

SliceOutputs = Union{
    HorizontalSlice,LongitudeDepthSlice,LatitudeDepthSlice,XDepthSlice,YDepthSlice
}

"""
$(TYPEDEF)

Free surface fields output.
    
$(TYPEDSIGNATURES)
    
## Details

Records two-dimensional free surface (height and barotropic velocity) fields.

$(TYPEDFIELDS)
"""
@kwdef struct FreeSurfaceFields{S} <: AbstractHorizontalModelOutput{S}
    "Time interval to record output at / s"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Depth averaged output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal fields corresponding to depth averaged model velocity and tracer
fields.

$(TYPEDFIELDS)
"""
@kwdef struct DepthAveraged{S} <: AbstractHorizontalModelOutput{S}
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

"""
$(TYPEDEF)

Meridional overturning circulation (MOC) stream function output.
    
$(TYPEDSIGNATURES)
    
## Details

Records latitude-depth or y-depth fields corresponding to stream function of meridional overturning
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

$(TYPEDFIELDS)
"""
@kwdef struct MOCStreamFunction{S} <: AbstractVerticalModelOutput{S}
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

"""
$(TYPEDEF)

Barotropic stream function output.
    
$(TYPEDSIGNATURES)
    
## Details

Records longitude-latitude or x-y fields corresponding to stream function of barotropic
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

$(TYPEDFIELDS)
"""
@kwdef struct BarotropicStreamFunction{S} <: AbstractHorizontalModelOutput{S}
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

LatitudeDepthOutputs = Union{LatitudeDepthSlice,MOCStreamFunction}
YDepthOutputs = Union{YDepthSlice,MOCStreamFunction}

@kwdef struct MOCStrength{S,T} <: AbstractScalarModelOutput{S}
    "y / m or latitude / ° coordinate to record output at"
    y_or_latitude::T
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(5day, window=5day)
end

@kwdef struct MeridionalHeatTransport{S,T} <: AbstractScalarModelOutput{S}
    "y / m or latitude / ° coordinate to record output at"
    y_or_latitude::T
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(5day, window=5day)
    "Sea water reference density / kg m⁻³"
    sea_water_density::T = 1026.0
    "Sea water specific heat capacity / J K⁻¹ kg⁻¹"
    sea_water_heat_capacity::T = 3991.0
end

@kwdef struct AverageKineticEnergy{S,M} <: AbstractScalarModelOutput{S}
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(5day, window=5day)
    "Mask defining region to average over or nothing if whole region should be used"
    region_mask::M = nothing
end

@kwdef struct HorizontallyAveragedTracers{S,M} <: AbstractModelOutput{S}
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(5day, window=5day)
    "Mask defining region to average over or nothing if whole region should be used"
    region_mask::M = nothing
end

ZeroOrOneDimensionalTimeSeriesOutputs = Union{
    HorizontallyAveragedTracers,AbstractScalarModelOutput
}

"""
Symbol label for output type to use in naming output file and registering output writer.

$(TYPEDSIGNATURES)
"""
function label(::AbstractModelOutput) end

label(output::HorizontalSlice) = Symbol("horizontal_slice_at_depth_$(output.depth)m")

function label(output::LongitudeDepthSlice)
    Symbol("longitude_depth_slice_at_latitude_$(output.latitude)deg")
end

function label(output::LatitudeDepthSlice)
    Symbol("latitude_depth_slice_at_longitude_$(output.longitude)deg")
end

label(output::XDepthSlice) = Symbol("x_depth_slice_at_y_$(output.y)m")
label(output::YDepthSlice) = Symbol("y_depth_slice_at_x_$(output.x)m")
label(::DepthAveraged) = :depth_averaged
label(::FreeSurfaceFields) = :free_surface_fields
label(::MOCStreamFunction) = :moc_stream_function
label(::BarotropicStreamFunction) = :barotropic_stream_function
label(output::MOCStrength) = Symbol("moc_strength_at_y_$(output.y_or_latitude)m_or_deg")
function label(output::MeridionalHeatTransport)
    Symbol("meridional_heat_transport_at_y_$(output.y_or_latitude)m_or_deg")
end

function label(base_label::Symbol, region_mask)
    isnothing(region_mask) ? base_label : Symbol("$(base_label)_$(summary(region_mask))")
end
label(output::AverageKineticEnergy) = label(:average_kinetic_energy, output.region_mask)
function label(output::HorizontallyAveragedTracers)
    label(:horizontally_averaged_tracers, output.region_mask)
end

"""
Spatial grid indices output type records fields at.

$(TYPEDSIGNATURES)
"""
indices(::AbstractModelOutput, grid::AbstractGrid) = (:, :, :)

indices(output::AbstractModelOutput, grid::ImmersedBoundaryGrid) = indices(output, grid.underlying_grid)

function indices(output::HorizontalSlice, grid::AbstractUnderlyingGrid)
    (:, :, clamp(searchsortedfirst(znodes(grid, Face()), output.depth), 1:grid.Nz))
end

function indices(output::LongitudeDepthSlice, grid::LatitudeLongitudeGrid)
    (:, clamp(searchsortedfirst(φnodes(grid, Face()), output.latitude), 1:grid.Ny), :)
end

function indices(output::LatitudeDepthSlice, grid::LatitudeLongitudeGrid)
    (clamp(searchsortedfirst(λnodes(grid, Face()), output.longitude), 1:grid.Nx), :, :)
end

function indices(output::XDepthSlice, grid::RectilinearGrid)
    (:, clamp(searchsortedfirst(ynodes(grid, Face()), output.y), 1:grid.Ny), :)
end

function indices(output::YDepthSlice, grid::RectilinearGrid)
    (clamp(searchsortedfirst(xnodes(grid, Face()), output.x), 1:grid.Nx), :, :)
end

function indices(
    output::Union{MOCStrength,MeridionalHeatTransport}, grid::LatitudeLongitudeGrid
)
    (:, clamp(searchsortedfirst(φnodes(grid, Face()), output.y_or_latitude), 1:grid.Ny), :)
end

function indices(output::Union{MOCStrength,MeridionalHeatTransport}, grid::RectilinearGrid)
    (:, clamp(searchsortedfirst(ynodes(grid, Face()), output.y_or_latitude), 1:grid.Ny), :)
end

"""
Named tuple of output variables (fields) to record for output type.
    
$(TYPEDSIGNATURES)
"""
function outputs(::AbstractModelOutput, model) end

outputs(::SliceOutputs, model) = merge(model.velocities, model.tracers)

function outputs(::DepthAveraged, model)
    NamedTuple(
        name => Field(Average(variable; dims=3)) for
        (name, variable) in pairs(merge(model.velocities, model.tracers))
    )
end

function outputs(
    ::FreeSurfaceFields,
    model::HydrostaticFreeSurfaceModel{
        TS,E,A,<:Union{ImplicitFreeSurface,ExplicitFreeSurface}
    },
) where {TS,E,A}
    (; η=model.free_surface.displacement)
end

function outputs(
    ::FreeSurfaceFields,
    model::HydrostaticFreeSurfaceModel{TS,E,A,<:SplitExplicitFreeSurface},
) where {TS,E,A}
    merge((; η=model.free_surface.displacement), model.free_surface.barotropic_velocities)
end

function outputs(::MOCStreamFunction, model)
    (;
        # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
        # units 10⁶ m³ s⁻¹ = Sv (Sverdrup)
        Ψᴹ=Field(
            CumulativeIntegral(
                Field(Integral(model.velocities.v * 1e-6; dims=1)); dims=3, reverse=true
            ),
        )
    )
end

function outputs(::BarotropicStreamFunction, model)
    (;
        # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
        # units 10⁶ m³ s⁻¹ = Sv (Sverdrup)
        Ψᴮ=Field(
            CumulativeIntegral(Field(Integral(model.velocities.v * 1e-6; dims=3)); dims=1)
        )
    )
end

function outputs(::MOCStrength, model)
    (;
        # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
        # units 10⁶ m³ s⁻¹ = Sv (Sverdrup)
        W=Field(
            Reduction(
                maximum!,
                Field(
                    CumulativeIntegral(
                        Field(Integral(model.velocities.v * 1e-6; dims=1));
                        dims=3,
                        reverse=true,
                    ),
                );
                dims=3,
            ),
        )
    )
end

function outputs(output::MeridionalHeatTransport, model)
    (;
        Q=Field(
            Integral(
                # Scale so integrated field is in units PW
                1e-15 *
                output.sea_water_density *
                output.sea_water_heat_capacity *
                model.velocities.v *
                model.tracers.T;
                dims=(1, 3),
            ),
        )
    )
end

function outputs(output::AverageKineticEnergy, model)
    (;
        eₖ=Field(
            Average(
                sum(v^2 / 2 for v in model.velocities);
                dims=(1, 2, 3),
                condition=output.region_mask,
            ),
        )
    )
end

function outputs(output::HorizontallyAveragedTracers, model)
    NamedTuple(
        name => Field(Average(tracer; dims=(1, 2), condition=output.region_mask)) for
        (name, tracer) in pairs(model.tracers)
    )
end

"""
Time schedule to record output type at.

$(TYPEDSIGNATURES)
"""
schedule(output::AbstractModelOutput) = output.schedule

"""
Filename to record outputs to.

$(TYPEDSIGNATURES)

## Details

For an output `output` a label computed using `label` function is appended on to
`stem`` and file extension is specified by `extension` added.
"""
function output_filename(stem::String, label::Symbol, extension::String)
    "$(stem)_$(label).$(extension)"
end

function output_filename(
    stem::String, output::AbstractModelOutput, extension::String="jld2"
)
    output_filename(stem, label(output), extension)
end

"""
    $(FUNCTIONNAME)(output, grid)

Label for horizontal axis for field heatmaps on grid `grid` for output `output`.
    
$(TYPEDSIGNATURES)
"""
function axis_xlabel end

axis_xlabel(::AbstractHorizontalModelOutput, ::LatitudeLongitudeGrid) = "Longitude λ / ᵒ"
axis_xlabel(::LongitudeDepthSlice, ::LatitudeLongitudeGrid) = "Longitude λ / ᵒ"
function axis_xlabel(::LatitudeDepthOutputs, ::LatitudeLongitudeGrid)
    "Latitude ϕ / ᵒ"
end
axis_xlabel(::AbstractHorizontalModelOutput, ::RectilinearGrid) = "x / m"
axis_xlabel(::XDepthSlice, ::RectilinearGrid) = "x / m"
axis_xlabel(::YDepthOutputs, ::RectilinearGrid) = "y / m"
function axis_xlabel(output::AbstractModelOutput, grid::ImmersedBoundaryGrid)
    (axis_xlabel(output, grid.underlying_grid))
end

function axis_xlabel(::ZeroOrOneDimensionalTimeSeriesOutputs, ::AbstractUnderlyingGrid)
    "Time / days"
end

"""
    $(FUNCTIONNAME)(output, grid)

Label for vertical axis for field heatmaps on grid `grid` for output `output`.
    
$(TYPEDSIGNATURES)
"""
function axis_ylabel end

axis_ylabel(::AbstractHorizontalModelOutput, ::LatitudeLongitudeGrid) = "Latitude ϕ / ᵒ"
axis_ylabel(::AbstractHorizontalModelOutput, ::RectilinearGrid) = "y / m"
function axis_ylabel(::AbstractVerticalModelOutput, ::AbstractUnderlyingGrid)
    "Depth z / m"
end
function axis_ylabel(output::AbstractModelOutput, grid::ImmersedBoundaryGrid)
    axis_ylabel(output, grid.underlying_grid)
end

axis_ylabel(::HorizontallyAveragedTracers, ::AbstractUnderlyingGrid) = "Depth z / m"

axis_ylabel(::AbstractScalarModelOutput, ::AbstractUnderlyingGrid) = ""

"""
Aspect ratio for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
axis_aspect_ratio(::AbstractModelOutput, ::AbstractGrid) = nothing
function axis_aspect_ratio(::AbstractHorizontalModelOutput, grid::LatitudeLongitudeGrid)
    AxisAspect(
        abs(-(extrema(λnodes(grid, Face()))...)) / abs(-(extrema(φnodes(grid, Face()))...))
    )
end
function axis_aspect_ratio(::AbstractHorizontalModelOutput, grid::RectilinearGrid)
    AxisAspect(
        abs(-(extrema(xnodes(grid, Face()))...)) / abs(-(extrema(ynodes(grid, Face()))...))
    )
end
function axis_aspect_ratio(output::AbstractModelOutput, grid::ImmersedBoundaryGrid)
    (axis_aspect_ratio(output, grid.underlying_grid))
end

"""
    $(FUNCTIONNAME)(output, grid)

Axis limits for field heatmaps on grid `grid` for output type `output`.
    
$(TYPEDSIGNATURES)
"""
function axis_limits end

function axis_xlimits(
    ::AbstractHorizontalModelOutput, grid::LatitudeLongitudeGrid, ::AbstractVector
)
    extrema(λnodes(grid, Face()))
end
function axis_xlimits(
    ::AbstractHorizontalModelOutput, grid::RectilinearGrid, ::AbstractVector
)
    extrema(xnodes(grid, Face()))
end
function axis_xlimits(::LongitudeDepthSlice, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(λnodes(grid, Face()))
end
function axis_xlimits(::LatitudeDepthOutputs, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(φnodes(grid, Face()))
end
function axis_xlimits(::XDepthSlice, grid::RectilinearGrid, ::AbstractVector)
    extrema(xnodes(grid, Face()))
end
function axis_xlimits(::YDepthSlice, grid::RectilinearGrid, ::AbstractVector)
    extrema(ynodes(grid, Face()))
end

function axis_xlimits(
    ::ZeroOrOneDimensionalTimeSeriesOutputs, ::AbstractUnderlyingGrid, times::AbstractVector
)
    extrema(times / 1day)
end

function axis_ylimits(
    ::AbstractHorizontalModelOutput, grid::LatitudeLongitudeGrid, ::AbstractVector
)
    extrema(φnodes(grid, Face()))
end
function axis_ylimits(
    ::AbstractHorizontalModelOutput, grid::RectilinearGrid, ::AbstractVector
)
    extrema(ynodes(grid, Face()))
end
function axis_ylimits(
    ::AbstractVerticalModelOutput, grid::AbstractUnderlyingGrid, ::AbstractVector
)
    extrema(znodes(grid, Face()))
end
function axis_ylimits(
    ::HorizontallyAveragedTracers, grid::AbstractUnderlyingGrid, ::AbstractVector
)
    extrema(znodes(grid, Face()))
end
function axis_ylimits(::AbstractScalarModelOutput, grid::AbstractUnderlyingGrid, ::AbstractVector)
    nothing
end
function axis_limits(
    output::AbstractModelOutput, grid::ImmersedBoundaryGrid, times::AbstractVector
)
    axis_limits(output, grid.underlying_grid, times)
end

function axis_limits(
    output::AbstractModelOutput, grid::AbstractUnderlyingGrid, times::AbstractVector
)
    (axis_xlimits(output, grid, times), axis_ylimits(output, grid, times))
end

@kwdef struct AutoVariableLimits{T}
    extrema_scale_factor::T = 0.8
    symmetrize::Bool = true
end

"""
Variable value limits for field color mapping.

$(TYPEDSIGNATURES)
"""
variable_limits(limits::Tuple{T,T}, ::FieldTimeSeries) where {T} = limits

function variable_limits(auto::AutoVariableLimits{T}, field_timeseries) where {T}
    limits = extrema(interior(field_timeseries)) .* auto.extrema_scale_factor
    auto.symmetrize ? (-maximum(abs.(limits)), maximum(abs.(limits))) : limits
end

struct VariablePlotConfiguration
    label::String
    unit::String
    color_map::Symbol
    limits::Union{Tuple,AutoVariableLimits}
end

const DEFAULT_VARIABLE_PLOT_CONFIGURATIONS = Dict{String,VariablePlotConfiguration}(
    "u" => VariablePlotConfiguration(
        "Zonal velocity u", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "v" => VariablePlotConfiguration(
        "Meridional velocity v", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "w" => VariablePlotConfiguration(
        "Vertical velocity w", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "T" => VariablePlotConfiguration("Temperature T", "°C", :thermal, (2.0, 30.0)),
    "S" => VariablePlotConfiguration("Salinity S", "g kg⁻¹", :haline, (32.0, 38.0)),
    "e" => VariablePlotConfiguration(
        "Turbulent kinetic energy e",
        "m² s⁻²",
        :amp,
        AutoVariableLimits(; symmetrize=false, extrema_scale_factor=0.5),
    ),
    "η" => VariablePlotConfiguration(
        "Free surface height η", "m", :balance, AutoVariableLimits()
    ),
    "U" => VariablePlotConfiguration(
        "Barotropic zonal velocity U", "m² s⁻¹", :balance, AutoVariableLimits()
    ),
    "V" => VariablePlotConfiguration(
        "Barotropic meridional velocity V", "m² s⁻¹", :balance, AutoVariableLimits()
    ),
    "Ψᴹ" => VariablePlotConfiguration(
        "MOC stream function Ψᴹ", "Sv", :balance, GyreInABox.AutoVariableLimits()
    ),
    "Ψᴮ" => VariablePlotConfiguration(
        "Barotropic stream function Ψᴮ", "Sv", :balance, GyreInABox.AutoVariableLimits()
    ),
    "W" => VariablePlotConfiguration(
        "MOC strength", "Sv", :balance, GyreInABox.AutoVariableLimits()
    ),
    "Q" => VariablePlotConfiguration(
        "Meridional heat transport", "PW", :balance, GyreInABox.AutoVariableLimits()
    ),
    "eₖ" => VariablePlotConfiguration(
        "Average kinetic energy",
        "m² s⁻²",
        :balance,
        GyreInABox.AutoVariableLimits(; symmetrize=false),
    ),
)

"""
Get unit associated with a variable.

$(SIGNATURES)
"""
unit(variable_name::String) = DEFAULT_VARIABLE_PLOT_CONFIGURATIONS[variable_name].unit
unit(variable::Symbol) = unit(string(variable))

"""
Register output writers for output types in `output_types` in `simulation` using
output filenames with stem `output_filename_stem`

$(SIGNATURES)
"""
function add_output_writers!(
    simulation::Oceananigans.Simulation, output_types::Tuple, output_filename_stem::String
)
    model = simulation.model
    # Create a grid on CPU if not already to avoid issues with computing output field
    # indices from grid using scalar operations on GPU grids
    grid = (
        isa(model.grid.architecture, CPU) ? model.grid : on_architecture(CPU(), model.grid)
    )

    for output in output_types
        simulation.output_writers[label(output)] = JLD2Writer(
            model,
            outputs(output, model);
            filename=output_filename(output_filename_stem, output),
            indices=indices(output, grid),
            schedule=schedule(output),
            overwrite_existing=true,
            with_halos=true,
        )
    end

    nothing
end

"""
Get properties for customizing plot axis rendering for `model_output` on `grid` and `times`.

$(SIGNATURES)
"""
function axis_properties(model_output, grid, times)
    (
        xlabel=axis_xlabel(model_output, grid),
        ylabel=axis_ylabel(model_output, grid),
        aspect=axis_aspect_ratio(model_output, grid),
        limits=axis_limits(model_output, grid, times),
    )
end

"""
Get ordered sequence of field variables to plot from `field_time_series`
based on those for which plot configurations are defined in
`variable_plot_configurations` and not excluded in `exclude_variables`.

$(SIGNATURES)
"""
function ordered_field_variables(
    field_timeseries, variable_plot_configurations, exclude_variables
)
    sort(
        tuple(
            setdiff(
                keys(field_timeseries) ∩ keys(variable_plot_configurations),
                exclude_variables,
            )...,
        );
        # Sort by tuple of variable name length and variable so that decorated variable
        # names appear after undecorated variables
        by=variable -> (length(variable), variable),
    )
end

"""
Get dicitionary of per-variable plot configurations by merging user-provided
overrides in `plot_configuration_overrides` with default configurations.

$(SIGNATURES)
"""
function get_variable_plot_configurations(plot_configuration_overrides)
    if isnothing(plot_configuration_overrides)
        DEFAULT_VARIABLE_PLOT_CONFIGURATIONS
    else
        merge(DEFAULT_VARIABLE_PLOT_CONFIGURATIONS, plot_configuration_overrides)
    end
end

"""
Compute dimensions of plot grid for `n_fields` fields with maximum number of
grid columns `max_columns`.

$(SIGNATURES)
"""
function plot_grid_dimensions(n_fields, max_columns)
    n_columns = min(n_fields, max_columns)
    n_rows = cld(n_fields, n_columns)
    return (n_columns, n_rows)
end

"""
Compute figure row and column indices for axis for field variable indexed by
`variable_index` for plot grid with `n_columns` columns starting at row
`row_offset`.

$(SIGNATURES)
"""
function plot_row_column_indices(variable_index, n_columns, row_offset)
    row = (variable_index - 1) ÷ n_columns + row_offset
    col = ((variable_index - 1) % n_columns) * 2 + 1
    return (row, col)
end

"""
Setup figure object of appropriate size for `n_rows` rows and `n_columns` of
axis objects each of size `(axis_width, axis_height)` plus a top margin of
`title_height` for inclusion of a figure title.

$(SIGNATURES)
"""
function setup_figure(n_rows, n_columns, axis_width, axis_height, title_height)
    Figure(;
        # CairoMakie defaults to px_per_unit=2 so manually adjust figure size here to
        # account for this - this is done in preference to changing px_per_unit using
        # CairoMakie.activate! to avoid persisting change after function exit
        size=(axis_width * n_columns / 2, (axis_height * n_rows + title_height) / 2),
        fontsize=12,
    )
end

abstract type AbstractPlotOutput end

"""
    $(FUNCTIONNAME)(axis, plot_output, field_timeseries, time_index, config)

Plots visual representation of `field_timeseries` appropriate for `plot_output`
output type on `axis`, optionally using observable `time_index` to index into
`field_timeseries` and with plot configuration options for variable represented
in `field_timeseries` specified in `config`.

$(TYPEDSIGNATURES)
"""
function plot_field_on_axis! end

"""
    $(FUNCTIONNAME)(plot_output, time_index, times)

Constructs figure title for `plot_output` output type optionally using information
about simulation times in `times` and observable time index `time_index`.

$(TYPEDSIGNATURES)
"""
function get_title end

"""
    $(FUNCTIONNAME)(
        plot_output,
        fig,
        times,
        time_index,
        output_filename_stem,
        model_output
    )

Saves plot file output for `plot_output` and `model_output` visualized on
Makie figure `fig` with model outputs recorded with `output_filename_step`, optionally
using simulation times `times` and observable `time_index`.

$(TYPEDSIGNATURES)
"""
function save_output end

"""
$(TYPEDEF)

Animated field plot output type.

## Details

Specifies recording an animation of model output fields recorded during a simulation.

$(TYPEDFIELDS)
"""
@kwdef struct AnimationPlotOutput <: AbstractPlotOutput
    "Frame rate (frames per second) to record animation at."
    frame_rate::Int = 10
    "Number of time indices to step through in field time series on each frame."
    frame_step::Int = 1
end

function plot_field_on_axis!(
    axis, ::AnimationPlotOutput, field_timeseries, time_index, config
)
    color_range = variable_limits(config.limits, field_timeseries)
    field = @lift field_timeseries[$time_index]
    heatmap!(axis, field; colormap=config.color_map, colorrange=color_range)
end

function get_title(::AnimationPlotOutput, time_index, times)
    @lift @sprintf("t = %s", prettytime(times[$time_index]))
end

function save_output(
    plot_output::AnimationPlotOutput,
    fig,
    times,
    time_index,
    output_filename_stem,
    model_output,
)
    frames = 1:plot_output.frame_step:length(times)

    output_file = output_filename(output_filename_stem, model_output, "mp4")
    @info "Recording an animation of $(label(model_output)) to $(output_file)..."

    CairoMakie.record(fig, output_file, frames; framerate=plot_output.frame_rate) do i
        (i % 10 == 0) && @printf "Plotting frame %i of %i\n" i frames[end]
        time_index[] = i
    end
end

"""
$(TYPEDEF)

Field temporal average plot output type.

## Details

Specifies plotting temporal average of model output fields recorded during a simulation.

The temporal averages are plotted as filled contour plots.

$(TYPEDFIELDS)
"""
@kwdef struct TemporalAveragePlotOutput{L<:Union{Int,AbstractVector}} <: AbstractPlotOutput
    """
    Either an integer specifying number of contour levels with range automatically
    determined or a vector of specific edge values to use.
    """
    levels::L = 10
end

function plot_field_on_axis!(
    axis, plot_output::TemporalAveragePlotOutput, field_timeseries, time_index, config
)
    # mean(field_timeseries, dims=4) errors on Oceananigans v105.2 so manually construct mean
    mean_field = sum(field_timeseries; dims=4) / length(field_timeseries)
    contourf!(axis, mean_field; colormap=config.color_map, levels=plot_output.levels)
end

function get_title(::TemporalAveragePlotOutput, time_index, times)
    @sprintf("Time average from\n%s to %s", prettytime(times[1]), prettytime(times[end]))
end

function save_output(
    ::TemporalAveragePlotOutput, fig, times, time_index, output_filename_stem, model_output
)
    save(output_filename(output_filename_stem, model_output, "svg"), fig)
end

"""
$(TYPEDEF)

Time series plot output type.

## Details

Specifies plotting time series of zero or one dimensional model output fields recorded during a simulation.

$(TYPEDFIELDS)
"""
@kwdef struct TimeSeriesPlotOutput{L<:Union{Int,AbstractVector}} <: AbstractPlotOutput
    """
    Either an integer specifying number of contour levels with range automatically
    determined or a vector of specific edge values to use (for 1D output fields).
    """
    levels::L = 10
end

function plot_field_on_axis!(
    axis, plot_output::TimeSeriesPlotOutput, field_timeseries, time_index, config
)
    times = field_timeseries.times / 1day
    dims = dim_x, dim_y, dim_z, dim_t = size(field_timeseries)
    indices = field_timeseries.indices
    indices = Tuple(
        s == 1 && indices[i] != Colon() ? first(indices[i]) : (s > 1 ? Colon() : 1)
        for (i, s) in enumerate(dims)
    )
    if dim_x == dim_y == dim_z == 1
        lines!(axis, times / 1day, field_timeseries[indices...])
        nothing
    else
        free_dimension_values= if dim_x > 1 && dim_y == dim_z == 1
            xnodes(field_timeseries)
        elseif dim_y > 1 && dim_x == dim_z == 1
            ynodes(field_timeseries)
        else
            znodes(field_timeseries)
        end
        contourf!(
            axis,
            times,
            free_dimension_values,
            field_timeseries[indices...]';
            colormap=config.color_map,
            levels=plot_output.levels,
        )
    end
end

get_title(::TimeSeriesPlotOutput, time_index, times) = nothing

function save_output(
    ::TimeSeriesPlotOutput, fig, times, time_index, output_filename_stem, model_output
)
    save(output_filename(output_filename_stem, model_output, "svg"), fig)
end

use_colorbar(::TimeSeriesPlotOutput) = true

"""
Create plot of fields recorded as model output.

$(SIGNATURES)

## Details

Generates and record to file with filename stem `output_filename_stem` an animation of
model outputs for output type `model_output` and for a model grid `grid`. The simulation
must have already been run for a model with this grid and with specified output type
active.

Animation is recorded with a frame rate of `frame_rate` frames per second, with fields
arranged on a grid with a maximum of `max_columns` columns, with the axis for each
field heatmap of size `(axis_width, axis_height)` in pixels and a further `title_height`
pixels allowed at the top of the figure for a title showing the simulation time. Each
frame of animation steps `frame_step` time indices through field time series.

By default the plot configurations for each field variable are taken from the
`GyreInABox.DEFAULT_VARIABLE_PLOT_CONFIGURATIONS` dictionary but these can be overridden
by passing a dictionary `plot_configuration_overrides` with keys corresponding to the
variable names to override. Specific variables to exclude from plot can be specified
in a tuple of variable names `exclude_variables`.
"""
function plot_output(
    plot_output_type::AbstractPlotOutput,
    output_filename_stem::String,
    model_output::AbstractModelOutput,
    grid::AbstractGrid;
    max_columns::Int=3,
    axis_width::Int=640,
    axis_height::Int=480,
    title_height::Int=40,
    exclude_variables::Tuple=(),
    plot_configuration_overrides::Union{Dict,Nothing}=nothing,
)
    filepath = output_filename(output_filename_stem, model_output)
    field_timeseries = FieldDataset(filepath).fields

    times = first(values(field_timeseries)).times

    axis_kwargs = axis_properties(model_output, grid, times)

    variable_plot_configurations = get_variable_plot_configurations(
        plot_configuration_overrides
    )

    field_variables = ordered_field_variables(
        field_timeseries, variable_plot_configurations, exclude_variables
    )

    n_columns, n_rows = plot_grid_dimensions(length(field_variables), max_columns)

    fig = setup_figure(n_rows, n_columns, axis_width, axis_height, title_height)

    time_index = Observable(1)

    title = get_title(plot_output_type, time_index, times)

    row_offset = isnothing(title) ? 1 : 2

    for (variable_index, variable) in enumerate(field_variables)
        config = variable_plot_configurations[variable]
        row, col = plot_row_column_indices(variable_index, n_columns, row_offset)
        axis = Axis(fig[row, col]; title="$(config.label) / $(config.unit)", axis_kwargs...)
        artist = plot_field_on_axis!(
            axis, plot_output_type, field_timeseries[variable], time_index, config
        )
        !isnothing(artist) && Colorbar(fig[row, col + 1], artist)
    end

    if !isnothing(title)
        fig[1, 1:end] = Label(fig, title; fontsize=20, tellwidth=false)
    end

    resize_to_layout!(fig)

    save_output(
        plot_output_type, fig, times, time_index, output_filename_stem, model_output
    )

    fig
end
