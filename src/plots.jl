"""
Get properties for customizing plot axis rendering for `variable` with 
`processor` on `grid` and `times`.

$(SIGNATURES)
"""
function axis_properties(variable::AbstractModelVariable, processor, grid, times)
    dimensions = spatial_dimensions(variable, processor)
    (
        xlabel=axis_xlabel(dimensions, grid),
        ylabel=axis_ylabel(dimensions, grid),
        aspect=axis_aspect_ratio(dimensions),
        limits=(
            axis_xlimits(dimensions, grid, times), axis_ylimits(dimensions, grid, times)
        ),
    )
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
        output_directory,
        output_filename_stem,
        model_output
    )

Saves plot file output for `plot_output` and `model_output` visualized on
Makie figure `fig` with model outputs recorded to a file in directory 
`output_directory` and with stem `output_filename_stem`, optionally
using simulation times `times` and observable `time_index`.

$(TYPEDSIGNATURES)
"""
function save_output end

"""
    $(FUNCTIONNAME)(plot_output, dimensions)

Indicates whether `plot_output` type is compatible with field with spatial dimensions
`dimensions` as a boolean.

$(TYPEDSIGNATURES)
"""
function is_compatible end

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
    axis,
    ::AnimationPlotOutput,
    field_timeseries,
    time_index,
    variable::AbstractModelVariable,
    limits::Tuple{<:Real,<:Real},
)
    field = @lift field_timeseries[$time_index]
    heatmap!(axis, field; colormap=color_map(variable), colorrange=limits)
end

function get_title(::AnimationPlotOutput, time_index, times)
    @lift @sprintf("t = %s", prettytime(times[$time_index]))
end

function save_output(
    plot_output::AnimationPlotOutput,
    fig,
    times,
    time_index,
    output_directory,
    output_filename_stem,
    model_output,
)
    frames = 1:plot_output.frame_step:length(times)

    output_file = joinpath(
        output_directory, output_filename(output_filename_stem, model_output, "mp4")
    )
    @info "Recording an animation of $(label(model_output)) to $(output_file)..."

    CairoMakie.record(fig, output_file, frames; framerate=plot_output.frame_rate) do i
        (i % 10 == 0) && @printf "Plotting frame %i of %i\n" i frames[end]
        time_index[] = i
    end
end

is_compatible(::AnimationPlotOutput, ::TwoSpatialDimensions) = true
is_compatible(::AnimationPlotOutput, ::ZeroOrOneSpatialDimensions) = false

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
    axis,
    plot_output::TemporalAveragePlotOutput,
    field_timeseries,
    time_index,
    variable::AbstractModelVariable,
    ::Tuple{<:Real,<:Real},
)
    # mean(field_timeseries, dims=4) errors on Oceananigans v105.2 so manually construct mean
    mean_field = sum(field_timeseries; dims=4) / length(field_timeseries)
    contourf!(axis, mean_field; colormap=color_map(variable), levels=plot_output.levels)
end

function get_title(::TemporalAveragePlotOutput, time_index, times)
    @sprintf("Time average from\n%s to %s", prettytime(times[1]), prettytime(times[end]))
end

function save_output(
    ::TemporalAveragePlotOutput,
    fig,
    times,
    time_index,
    output_directory,
    output_filename_stem,
    model_output,
)
    save(
        joinpath(
            output_directory, output_filename(output_filename_stem, model_output, "svg")
        ),
        fig,
    )
end

is_compatible(::TemporalAveragePlotOutput, ::TwoSpatialDimensions) = true
is_compatible(::TemporalAveragePlotOutput, ::ZeroOrOneSpatialDimensions) = false

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
    axis,
    plot_output::TimeSeriesPlotOutput,
    field_timeseries,
    time_index,
    variable::AbstractModelVariable,
    ::Tuple{<:Real,<:Real},
)
    times = field_timeseries.times / 1day
    spatial_dims..., _ = dim_x, dim_y, dim_z, dim_t = size(field_timeseries)
    indices = field_timeseries.indices
    indices = Tuple(
        s == 1 && indices[i] != Colon() ? first(indices[i]) : (s > 1 ? Colon() : 1) for
        (i, s) in enumerate(spatial_dims)
    )
    if dim_x == dim_y == dim_z == 1
        lines!(axis, times, [field_timeseries[i][indices...] for i in 1:length(times)])
        nothing
    else
        free_dimension_values = if dim_x > 1 && dim_y == dim_z == 1
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
            hcat((interior(field_timeseries[i])[indices...] for i in 1:length(times))...)';
            colormap=color_map(variable),
            levels=plot_output.levels,
        )
    end
end

get_title(::TimeSeriesPlotOutput, time_index, times) = nothing

function save_output(
    ::TimeSeriesPlotOutput,
    fig,
    times,
    time_index,
    output_directory,
    output_filename_stem,
    model_output,
)
    save(
        joinpath(
            output_directory, output_filename(output_filename_stem, model_output, "svg")
        ),
        fig,
    )
end

is_compatible(::TimeSeriesPlotOutput, ::TwoSpatialDimensions) = false
is_compatible(::TimeSeriesPlotOutput, ::ZeroOrOneSpatialDimensions) = true

"""
$(SIGNATURES)

Create plot of fields recorded as model output.

## Details

Generates and record to file in `output_directory` with filename stem `output_filename_stem`
a plot of type `plot_output_type` for model outputs for output type `model_output` and for
a model grid `grid`. The simulation must have already been run for a model with this grid
and with specified output type active and recorded with extension `output_file_extension`
specific to the output writer type.

Fields are arranged on a grid with a maximum of `max_columns` columns, with the axis for
each field heatmap of size `(axis_width, axis_height)` in pixels and a further `title_height`
pixels allowed at the top of the figure for a title.

By default the value limits for each field variable are automatically computed from
the extrema of the field data but these can be overridden by passing a named tuple 
`variable_limits` with keys corresponding to the variable (short) names to override. 
Specific variables to exclude from plot can be specified in a tuple of variable names
`exclude_variables`.

The backend used for loading data from file can be set using `backend` (defaults to
loading all data in memory). A subset of times at which output was recorded at
which should be used to produce plot can be specified using `times` - the default
is to use all recorded times.
"""
function plot_output(
    plot_output_type::AbstractPlotOutput,
    output_directory::String,
    output_filename_stem::String,
    model_output::ModelOutput,
    grid::AbstractGrid;
    output_file_extension::String="jld2",
    max_columns::Int=3,
    axis_width::Int=640,
    axis_height::Int=480,
    title_height::Int=40,
    exclude_variables::Tuple=(),
    variable_limits::Union{NamedTuple,Nothing}=nothing,
    backend::Union{InMemory,OnDisk}=InMemory(),
    times::Union{AbstractVector,Nothing}=nothing,
)
    filepath = joinpath(
        output_directory,
        output_filename(output_filename_stem, model_output, output_file_extension),
    )

    variables = filter(v -> short_name(variable) ∉ exclude_variables, output.variables)

    field_data = Dict{String,FieldTimeSeries}(
        short_name(variable) =>
            FieldTimeSeries(filepath, short_name(variable); backend, times) for
        variable in variable
    )

    times = first(values(field_data)).times

    n_columns, n_rows = plot_grid_dimensions(length(variables), max_columns)

    fig = setup_figure(n_rows, n_columns, axis_width, axis_height, title_height)

    time_index = Observable(1)

    title = get_title(plot_output_type, time_index, times)

    row_offset = isnothing(title) ? 1 : 2

    for (variable_index, variable) in enumerate(variables)
        axis_kwargs = axis_properties(variable, model_output.processor, grid, times)
        row, col = plot_row_column_indices(variable_index, n_columns, row_offset)
        axis = Axis(
            fig[row, col];
            title="$(long_name(variable)) / $(units(variable))",
            axis_kwargs...,
        )
        name = short_name(variable)
        limits = get(variable_limits, Symbol(name), nothing)
        if isnothing(limits)
            limits = auto_variable_limits(variable, field_data[name])
        end
        artist = plot_field_on_axis!(
            axis, plot_output_type, field_data[name], time_index, variable, limits
        )
        !isnothing(artist) && Colorbar(fig[row, col + 1], artist)
    end

    if !isnothing(title)
        fig[1, 1:end] = Label(fig, title; fontsize=20, tellwidth=false)
    end

    resize_to_layout!(fig)

    save_output(
        plot_output_type,
        fig,
        times,
        time_index,
        output_directory,
        output_filename_stem,
        model_output,
    )

    fig
end