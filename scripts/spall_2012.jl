using ArgParse
using GyreInABox
using JLD2
using Oceananigans
using Oceananigans.Units
using CUDA
using Dates: now, format

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--surface-temperature-restoring-strength", "-R"
            help = "Surface temperature restoring strength / W m⁻² K⁻¹"
            arg_type = Float64
            default = 20.
        "--northern-basin-surface-evaporation", "-E"
            help = "Surface net evaporation-precipitation in northern basin above sill / m s⁻¹"
            arg_type = Float64
            default = -2e-8
        "--sill-height", "-H"
            help = "Sill height / m"
            arg_type = Float64
            default = 1000.
        "--simulation-years", "-Y"
            help = "Number of simulated year to run for"
            arg_type = Float64
            default = 100.
        "--output-interval-days", "-I"
            help = "Interval at which to record outputs at in simulated days"
            arg_type = Float64
            default = 30.
        "--grid-size", "-G"
            help = "Grid dimensions in x, y and depth"
            nargs = 3
            arg_type = Int
            default = [200, 400, 30]
        "--cpu"
            help = "Run on CPU (rather than GPU, the default)"
            action = :store_true
        "--use-eddy-closure"
            help = "Use a dynamic Smagorinsky eddy closure"
            action = :store_true
    end

    return parse_args(s)
end

function main()

    args = parse_commandline()

    @info "Parsed args:"
    for (arg, val) in args
        @info "  $arg: $val"
    end

    architecture = args["cpu"] ? CPU() : GPU()

    output_interval = args["output-interval-days"] * 1day

    parameters = Spall2011Parameters(;
        grid_size=Tuple(args["grid-size"]),
        surface_temperature_restoring_strength=args["surface-temperature-restoring-strength"],
        northern_basin_surface_evaporation=args["northern-basin-surface-evaporation"],
        sill_height=args["sill-height"],
        use_eddy_closure=args["use-eddy-closure"],
    )

    mask = GyreInABox.northern_basin_mask(parameters)

    configuration = SimulationConfiguration(
        architecture=architecture,
        simulation_time=args["simulation-years"]*365day,
        initial_timestep=10minute,
        maximum_timestep=60minute,
        wizard_max_change=1.1,
        wizard_update_interval=10,
        output_filename="spall_2012_gyre_model",
        output_types=(
            HorizontalSlice(schedule=TimeInterval(output_interval)),
            HorizontalSlice(depth=parameters.bottom_depth + parameters.sill_height, schedule=TimeInterval(output_interval)),
            XDepthSlice(y_or_latitude=parameters.sill_center_y, schedule=TimeInterval(output_interval)),
            YDepthSlice(x_or_longitude=parameters.domain_size_x / 2, schedule=TimeInterval(output_interval)),
            FreeSurfaceFields(schedule=TimeInterval(output_interval)),
            MOCStreamFunction(),
            BarotropicStreamFunction(),
            MOCStrength(y_or_latitude=parameters.sill_center_y),
            MeridionalHeatTransport(y_or_latitude=parameters.sill_center_y),
            AverageKineticEnergy(region_mask=mask),
            HorizontallyAveragedTracers(region_mask=mask),
        ),
        progress_message_interval=1000
    )
    
    timestamp = format(now(), "yyyy-mm-dd_HH-MM-SS")

    run_directory = mkdir("$timestamp-run")
    cd(run_directory)

    jldsave("parameters.jld2"; parameters)
    jldsave("configuration.jld2"; configuration)

    run_simulation(parameters, configuration)

end

main()
