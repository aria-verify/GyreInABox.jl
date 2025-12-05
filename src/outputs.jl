abstract type AbstractOutput{S} end

abstract type AbstractHorizontalOutput{S} <: AbstractOutput{S} end

abstract type AbstractVerticalOutput{S} <: AbstractOutput{S} end

abstract type AbstractLatitudeDepthOutput{S} <: AbstractVerticalOutput{S} end

"""
$(TYPEDEF)

Horizontal (latitude-longitude) slice output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal slices through model velocity and tracer fields at specified depth.

$(TYPEDFIELDS)
"""
@kwdef struct HorizontalSlice{S, T} <: AbstractHorizontalOutput{S}
    "Depth of slice / m"
    depth::T = 0.
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
@kwdef struct LongitudeDepthSlice{S, T} <: AbstractVerticalOutput{S}
    "Latitude of slice / °"
    latitude::T = 45.
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
@kwdef struct LatitudeDepthSlice{S, T} <: AbstractLatitudeDepthOutput{S}
    "Longitude of slice / °"
    longitude::T = 30.
    "Schedule to record output at"
    schedule::S = TimeInterval(1day)
end

SliceOutputType = Union{HorizontalSlice, LongitudeDepthSlice, LatitudeDepthSlice}

"""
$(TYPEDEF)

Free surface fields output.
    
$(TYPEDSIGNATURES)
    
## Details

Records two-dimensional free surface (height and barotropic velocity) fields.

$(TYPEDFIELDS)
"""
@kwdef struct FreeSurfaceFields{S} <: AbstractHorizontalOutput{S}
    "Time interval to record output at / s"
    schedule::S = TimeInterval(1day)
end

"""
$(TYPEDEF)

Depth and time averaged output.
    
$(TYPEDSIGNATURES)
    
## Details

Records horizontal fields corresponding to depth averaged model velocity and tracer
fields.

$(TYPEDFIELDS)
"""
@kwdef struct DepthAveraged{S} <: AbstractHorizontalOutput{S}
    "Schedule to record output at"
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

"""
$(TYPEDEF)

Meridional overturning circulation (MOC) stream function output.
    
$(TYPEDSIGNATURES)
    
## Details

Records latitude-depth fields corresponding to stream function of meridional overturning
circulation - computed here as vertically accumulated - that is cumulative vertical
integral with respect to depth - of zonally integrated meridional velocity component:

```math
\\Psi^M(\\varphi, z, t) = 
\\int_{0}^z \\int_{\\lambda_W}^{\\lambda_E} 
  v(\\lambda, \\varphi, z', t)
\\,\\mathrm{d}\\lambda \\,\\mathrm{d}z'
```

The outputted field is scaled to be in sverdrup (10⁶ m³ s⁻¹) units.

$(TYPEDFIELDS)
"""
@kwdef struct MOCStreamFunction{S} <: AbstractLatitudeDepthOutput{S}
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

"""
$(TYPEDEF)

Barotropic stream function output.
    
$(TYPEDSIGNATURES)
    
## Details

Records longitude-latitude fields corresponding to stream function of barotropic
velocity - computed here as zonally accumulated - that is cumulative integral with
respect to longitude - of depth integrated meridional velocity component:

```math
\\Psi^B(\\lambda, \\varphi, t) = 
\\int_{\\lambda_W}^{\\lambda}  \\int_{z_B}^{z_S} 
  v(\\lambda', \\varphi, z, t)
\\,\\mathrm{d}z \\,\\mathrm{d}\\lambda'
```

The outputted field is scaled to be in sverdrup (10⁶ m³ s⁻¹) units.

$(TYPEDFIELDS)
"""
@kwdef struct BarotropicStreamFunction{S} <: AbstractHorizontalOutput{S}
    schedule::S = AveragedTimeInterval(30day, window=30day)
end

"""
Symbol label for output type to use in naming output file and registering output writer.

$(TYPEDSIGNATURES)
"""
function label(::AbstractOutput) end

label(output::HorizontalSlice) = Symbol("horizontal_slice_at_depth_$(output.depth)m")
label(output::LongitudeDepthSlice) = Symbol(
    "longitude_depth_slice_at_latitude_$(output.latitude)deg"
)
label(output::LatitudeDepthSlice) = Symbol(
    "latitude_depth_slice_at_longitude_$(output.longitude)deg"
)
label(::DepthAveraged) = :depth_averaged
label(::FreeSurfaceFields) = :free_surface_fields
label(::MOCStreamFunction) = :moc_stream_function
label(::BarotropicStreamFunction) = :barotropic_stream_function

"""
Spatial grid indices output type records fields at.

$(TYPEDSIGNATURES)
"""
indices(::AbstractOutput, grid) = (:, :, :)

indices(output::HorizontalSlice, grid) = (
    :, :, clamp(searchsortedfirst(znodes(grid, Face()), output.depth), 1:grid.Nz)
)
indices(output::LongitudeDepthSlice, grid) = (
    :, clamp(searchsortedfirst(φnodes(grid, Face()), output.latitude), 1:grid.Ny), :
)
indices(output::LatitudeDepthSlice, grid) = (
    clamp(searchsortedfirst(λnodes(grid, Face()), output.longitude), 1:grid.Nx), :, :
)

"""
Named tuple of output variables (fields) to record for output type.
    
$(TYPEDSIGNATURES)
"""
function outputs(::AbstractOutput, model) end

outputs(::SliceOutputType, model) = merge(model.velocities, model.tracers)
outputs(::DepthAveraged, model) = NamedTuple(
    name => Field(Average(variable, dims=3))
    for (name, variable) in pairs(merge(model.velocities, model.tracers))
)
outputs(::FreeSurfaceFields, model) = merge(
    (; η=model.free_surface.η),
    model.free_surface.barotropic_velocities
)
outputs(::MOCStreamFunction, model) = (;
    # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
    # units 10⁶ m³ s⁻¹ = Sv (Sverdrup)
    Ψᴹ = Field(
        CumulativeIntegral(

            Field(Integral(model.velocities.v * 1e-6, dims=1)), dims=3, reverse=true
        )
    )
)
outputs(::BarotropicStreamFunction, model) = (;
    # Scale velocities by 1 / 10⁶ so doubly spatially integrated field is in
    # units 10⁶ m³ s⁻¹ = Sv (Sverdrup)
    Ψᴮ = Field(
        CumulativeIntegral(
            Field(Integral(model.velocities.v * 1e-6, dims=3)), dims=1
        )
    )
)

"""
Time schedule to record output type at.

$(TYPEDSIGNATURES)
"""
schedule(output::AbstractOutput) = output.schedule

"""
Filename to record outputs to.

$(TYPEDSIGNATURES)

## Details

For an output type `output` a label computed using `label` function is appended on to
`stem`` and file extension is specified by `extension` added.
"""
output_filename(
    stem::String, label::Symbol, extension::String
) = "$(stem)_$(label).$(extension)"
output_filename(
    stem::String,
    output::AbstractOutput,
    extension::String="jld2"
) = output_filename(stem, label(output), extension)

"""
Label for horizontal axis for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_xlabel(::AbstractOutput) end

axis_xlabel(::AbstractHorizontalOutput) = "Longitude λ / ᵒ"
axis_xlabel(::LongitudeDepthSlice) = "Longitude λ / ᵒ"
axis_xlabel(::AbstractLatitudeDepthOutput) = "Latitude ϕ / ᵒ"

"""
Label for vertical axis for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_ylabel(::AbstractOutput) end

axis_ylabel(::AbstractHorizontalOutput) = "Latitude ϕ / ᵒ"
axis_ylabel(::AbstractVerticalOutput) = "Depth z / m"


"""
Aspect ratio for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
axis_aspect_ratio(::AbstractGrid, ::AbstractOutput) = AxisAspect(1)
axis_aspect_ratio(grid::AbstractGrid, ::AbstractHorizontalOutput) = AxisAspect(
    abs(-(extrema(λnodes(grid, Face()))...)) / abs(-(extrema(φnodes(grid, Face()))...))
)

"""
Axis limits for field heatmaps.
    
$(TYPEDSIGNATURES)
"""
function axis_limits(::AbstractGrid, ::AbstractOutput) end

axis_limits(grid::AbstractGrid, ::AbstractHorizontalOutput) = (
    extrema(λnodes(grid, Face())), extrema(φnodes(grid, Face()))
)
axis_limits(grid::AbstractGrid, ::LongitudeDepthSlice) = (
    extrema(λnodes(grid, Face())), extrema(znodes(grid, Face()))
)
axis_limits(grid::AbstractGrid, ::AbstractLatitudeDepthOutput) = (
    extrema(φnodes(grid, Face())), extrema(znodes(grid, Face()))
)

@kwdef struct AutoVariableLimits{T}
    extrema_scale_factor::T = 0.8
    symmetrize::Bool = true
end

"""
Variable value limits for field color mapping.

$(TYPEDSIGNATURES)
"""
variable_limits(limits::Tuple{T, T}, ::FieldTimeSeries) where {T} = limits

function variable_limits(auto::AutoVariableLimits{T}, field_timeseries) where {T}
    limits = extrema(interior(field_timeseries)) .* auto.extrema_scale_factor
    auto.symmetrize ? (-maximum(abs.(limits)), maximum(abs.(limits))) : limits 
end

struct VariablePlotConfiguration
    label::String
    unit::String
    color_map::Symbol
    limits::Union{Tuple, AutoVariableLimits}
end

const DEFAULT_VARIABLE_PLOT_CONFIGURATIONS = Dict{String, VariablePlotConfiguration}(
    "u" => VariablePlotConfiguration(
        "Zonal velocity u", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "v" => VariablePlotConfiguration(
        "Meridional velocity v", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "w" => VariablePlotConfiguration(
        "Vertical velocity w", "m s⁻¹", :balance, AutoVariableLimits()
    ),
    "T" => VariablePlotConfiguration(
        "Temperature T", "°C", :thermal, (2., 30.)
    ),
    "S" => VariablePlotConfiguration(
        "Salinity S", "g kg⁻¹", :haline, (32., 38.)
    ),
    "e" => VariablePlotConfiguration(
        "Turbulent kinetic energy e", 
        "m² s⁻²",
        :amp, 
        AutoVariableLimits(symmetrize=false, extrema_scale_factor=0.5)
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
)

"""
Get unit associated with a variable.

$(SIGNATURES)
"""
unit(variable_name::String) = DEFAULT_VARIABLE_PLOT_CONFIGURATIONS[variable_name].unit
unit(variable::Symbol) = unit(string(variable))

"""
Register output writers for output types specified in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_output_writers!(
    simulation::Oceananigans.Simulation,
    output_types::Tuple,
    output_filename_stem::String
)
    model = simulation.model
    # Create a grid on CPU if not already to avoid issues with computing output field
    # indices from grid using scalar operations on GPU grids
    grid = (
        isa(model.grid.architecture, CPU) 
        ? model.grid 
        : on_architecture(CPU(), model.grid)
    )

    for output_type in output_types
        simulation.output_writers[label(output_type)] = JLD2Writer(
            model,
            outputs(output_type, model),
            filename=output_filename(output_filename_stem, output_type),
            indices=indices(output_type, grid),
            schedule=schedule(output_type),
            overwrite_existing=true,
            with_halos=true,
        )
    end

    nothing
end


"""
Record animation of fields recorded as model output.

$(SIGNATURES)

## Details

Generates and record to file with filename stem `output_filename_stem` an animation of
model outputs for output type `output_type` and for a model grid `grid`. The simulation
must have already been run for a model with this grid and with specified output type
active.

Animation is recorded with a frame rate of `frame_rate` frames per second, with fields
arranged on a grid with a maximum of `max_columns` columns, with the axis for each
field heatmap of size `(axis_width, axis_height)` in pixels and a further `title_height`
pixels allowed at the top of the figure for a title showing the simulation time.

By default the plot configurations for each field variable are taken from the
`GyreInABox.DEFAULT_VARIABLE_PLOT_CONFIGURATIONS` dictionary but these can be overridden
by passing a dictionary `plot_configuration_overrides` with keys corresponding to the
variable names to override. Specific variables to exclude from plot can be specified
in a tuple of variable names `exclude_variables`.
"""
function record_animation(
    output_filename_stem::String,
    output_type::AbstractOutput,
    grid::AbstractGrid;
    frame_rate::Int = 10,
    max_columns::Int = 3,
    axis_width::Int = 640,
    axis_height::Int = 480,
    title_height::Int = 40,
    exclude_variables::Tuple = (),
    plot_configuration_overrides::Union{Dict, Nothing} = nothing
)
    filepath = output_filename(output_filename_stem, output_type)
    field_timeseries = FieldDataset(filepath).fields
    
    axis_kwargs = (
        xlabel=axis_xlabel(output_type),
        ylabel=axis_ylabel(output_type),
        aspect=axis_aspect_ratio(grid, output_type),
        limits=axis_limits(grid, output_type)
    )
    
    variable_plot_configurations = (
        isnothing(plot_configuration_overrides) 
        ? DEFAULT_VARIABLE_PLOT_CONFIGURATIONS 
        : merge(DEFAULT_VARIABLE_PLOT_CONFIGURATIONS, plot_configuration_overrides)
    )
    
    field_variables = sort(
        tuple(
            setdiff(
                keys(field_timeseries) ∩ keys(variable_plot_configurations),
                exclude_variables
            )...
        ),
        # Sort by tuple of variable name length and variable so that decorated variable
        # names appear after undecorated variables
        by=variable -> (length(variable), variable)
    )
    
    n_fields = length(field_variables)
    n_columns = min(n_fields, max_columns)
    n_rows = cld(n_fields, n_columns)
    
    fig = Figure(
        # CairoMakie defaults to px_per_unit=2 so manually adjust figure size here to
        # account for this - this is done in preference to changing px_per_unit using
        # CairoMakie.activate! to avoid persisting change after function exit
        size=(axis_width * n_columns / 2, (axis_height * n_rows + title_height) / 2),
        fontsize=12,
    )
    
    times = first(values(field_timeseries)).times
    time_index = Observable(1)
    
    for (variable_index, variable) in enumerate(field_variables)
        config = variable_plot_configurations[variable]
        color_range = variable_limits(config.limits, field_timeseries[variable])
        field = @lift field_timeseries[variable][$time_index]
        row = (variable_index - 1) ÷ n_columns + 2
        col = ((variable_index - 1) % n_columns) * 2 + 1
        axis = Axis(fig[row, col]; title=config.label, axis_kwargs...)
        hmap = heatmap!(axis, field; colormap=config.color_map, colorrange=color_range)
        Colorbar(fig[row, col + 1], hmap; label=config.unit)
    end
    
    title = @lift @sprintf("t = %s", prettytime(times[$time_index]))
    
    fig[1, 1:n_columns * 2] = Label(fig, title, fontsize=20, tellwidth=false)

    frames = 1:length(times)

    output_file = output_filename(output_filename_stem, output_type, "mp4")
    @info "Recording an animation of $(label(output_type)) to $(output_file)..."

    CairoMakie.record(fig, output_file, frames, framerate=frame_rate) do i
        (i % 10 == 0) && @printf "Plotting frame %i of %i\n" i frames[end]
        time_index[] = i
    end

    fig
end
