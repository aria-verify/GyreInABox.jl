using ArgParse
using GyreInABox
using JLD2
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

    outputs = Dict()

    for output_type in (
        HorizontallyAveragedTracers,
        AverageKineticEnergy,
        MOCStrength,
        MeridionalHeatTransport,
    )
        output_type_index = findfirst(ot -> ot isa output_type, configuration.output_types)
        outputs[output_type] = if !isnothing(output_type_index)
            configuration.output_types[output_type_index]
        else
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

    if !isnothing(outputs[HorizontallyAveragedTracers])
        averaged_tracers_output = outputs[HorizontallyAveragedTracers]
        S_timeseries = FieldTimeSeries(
            "$(args["run-output-path"])/$(GyreInABox.output_filename(configuration.output_filename, averaged_tracers_output))",
            "S";
            backend=InMemory(),
            times=plot_times,
        )
        T_timeseries = FieldTimeSeries(
            "$(args["run-output-path"])/$(GyreInABox.output_filename(configuration.output_filename, averaged_tracers_output))",
            "T";
            backend=InMemory(),
            times=plot_times,
        )
        ρ_timeseries = (
            parameters.sea_water_density .+ stack(
                (parameters.haline_contraction_coefficient * S_timeseries[t] - parameters.thermal_expansion_coefficient * T_timeseries[t])[
                    1, 1, 1:parameters.grid_size[3]
                ] for t in 1:length(S_timeseries)
            )
        )
        times = S_timeseries.times / 365days
        depths = znodes(S_timeseries)

        col = 1
        for (label, timeseries, cmap) in [
            (
                "Salinity / g kg⁻¹",
                S_timeseries[1, 1, 1:parameters.grid_size[3], 1:length(times)],
                :haline,
            ),
            (
                "Temperature / °C",
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

    for (output_type, field_name, label) in (
        (AverageKineticEnergy, "eₖ", "Average kinetic energy / m²s⁻²"),
        (MOCStrength, "W", "MOC strength / Sv"),
        (MeridionalHeatTransport, "Q", "Meridional heat transport / TW"),
    )
        if !isnothing(outputs[output_type])
            output = outputs[output_type]
            timeseries = FieldTimeSeries(
                "$(args["run-output-path"])/$(GyreInABox.output_filename(configuration.output_filename, output))",
                field_name;
                backend=InMemory(),
                times=plot_times,
            )
            times = timeseries.times / 365days
            ax = Axis(
                fig[row, col:(col + 1)];
                title=label,
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

    output_path = "$(args["run-output-path"])/$(configuration.output_filename)_$(args["output-filename-suffix"]).svg"

    @info "Writing output figure to $(output_path)"
    save(output_path, fig)
end

main()