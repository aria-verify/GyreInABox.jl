using ArgParse
using GyreInABox
using JLD2
using NCDatasets
using Zarr
using Oceananigans
using Oceananigans.Units
using CairoMakie

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "run-output-path"
            help = "Path to directory containing run outputs"
            arg_type = String
            required = true
        "--plot-time-range", "-t"
            help = "Time range for plots in format (start, step, end) in days"
            arg_type = Float64
            nargs = 3
        "--output-filename-suffix", "-o"
            help = "Output filename suffix for summary plot"
            arg_type = String
            default = "summary"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    parameters = load("$(args["run-output-path"])/parameters.jld2")["parameters"]
    configuration = load("$(args["run-output-path"])/configuration.jld2")["configuration"]

    outputs = Dict{GyreInABox.AbstractModelVariable, GyreInABox.ModelOutput}()

    variables = (
        S=Salinity(),
        T=Temperature(),
        eₖ=SpecificKineticEnergy(),
        Ψᴹ=MOCStreamFunction(),
        Q=NorthwardHeatTransport(),
    )

    for variable in values(variables)
        output_type_index = findfirst(
            ot -> (
                variable in ot.variables && 
                GyreInABox.spatial_dimensions(
                    variable, ot.processor
                ) isa GyreInABox.ZeroOrOneSpatialDimensions
            ), 
            configuration.output_types
        )
        outputs[variable] = if !isnothing(output_type_index)
            configuration.output_types[output_type_index]
        else
            @info "Model output for variable $(variable) not found in configuration"
            nothing
        end
    end

    plot_times = if !isempty(args["plot-time-range"])
        start, step, stop = args["plot-time-range"] .* 1day
        range(start, stop; step)
    else
        nothing
    end

    fig = Figure(; size=(1200, 800))

    row = 2

    function load_field_time_series(variable)
        FieldTimeSeries(
                joinpath(
                    args["run-output-path"],
                    GyreInABox.output_filename(
                        configuration.output_filename_stem, 
                        outputs[variable],
                        GyreInABox.extension(configuration.output_writer_type)
                    )
                ),
                GyreInABox.short_name(variable);
                backend=InMemory(),
                times=plot_times,
        )
    end

    if !isnothing(outputs[variables.S]) && !isnothing(outputs[variables.T])
        S_timeseries = load_field_time_series(variables.S)
        T_timeseries = load_field_time_series(variables.T)
        ρ_timeseries = (
            parameters.sea_water_density .+ stack(
                (parameters.haline_contraction_coefficient * S_timeseries[t] - parameters.thermal_expansion_coefficient * T_timeseries[t])[
                    1, 1, 1:parameters.grid_size[3]
                ] for t in 1:length(S_timeseries)
            )
        )
        times = collect(S_timeseries.times / 365days)
        depths = znodes(S_timeseries)

        col = 1
        for (label, timeseries, cmap) in [
            (
                "$(GyreInABox.long_name(variables.S)) / $(GyreInABox.units(variables.S))",
                S_timeseries[1, 1, 1:parameters.grid_size[3], 1:length(times)],
                :haline,
            ),
            (
                "$(GyreInABox.long_name(variables.T)) / $(GyreInABox.units(variables.T))",
                T_timeseries[1, 1, 1:parameters.grid_size[3], 1:length(times)],
                :thermal,
            ),
            ("Density / kg m⁻³", ρ_timeseries, :dense),
        ]
            ax = Axis(
                fig[row, col];
                title=label,
                xlabel="Time / years",
                ylabel="Depth / m",
                limits=(extrema(times), extrema(depths)),
            );
            col += 1
            cf = contourf!(ax, times, depths, timeseries'; colormap=cmap)
            Colorbar(fig[row, col], cf)
            col += 1
        end

        row += 1
    end

    col = 1

    for variable in (variables.eₖ, variables.Ψᴹ, variables.Q)
        if !isnothing(outputs[variable])
            timeseries = load_field_time_series(variable)
            times = collect(timeseries.times / 365days)
            ax = Axis(
                fig[row, col:(col + 1)];
                title="$(GyreInABox.long_name(variable)) / $(GyreInABox.units(variable))",
                xlabel="Time / years",
                limits=(extrema(times), nothing),
            )
            dims = size(timeseries)
            indices = timeseries.indices
            indices = Tuple(
                s == 1 && indices[i] != Colon() ? first(indices[i]) : (s > 1 ? Colon() : 1)
                for (i, s) in enumerate(dims)
            )
            lines!(ax, times, timeseries[indices...])
            col += 2
        end
    end

    fig[1, 1:end] = Label(
        fig,
        "Γ = $(parameters.surface_temperature_restoring_strength) Wm⁻²K⁻¹ and " *
        "E = $(parameters.northern_basin_surface_evaporation / 1e-8)×10⁻⁸ ms⁻¹";
        fontsize=20,
        tellwidth=false,
    )

    resize_to_layout!(fig)

    output_path = joinpath(
        args["run-output-path"],
        GyreInABox.output_filename(
            configuration.output_filename_stem,
            args["output-filename-suffix"],
            "svg"
        )
    )

    @info "Writing output figure to $(output_path)"
    save(output_path, fig)
end

main()