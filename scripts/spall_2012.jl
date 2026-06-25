using ArgParse
using GyreInABox
using JLD2
using Oceananigans
using Oceananigans.DistributedComputations
using Oceananigans.Units
using CUDA
using MPI
using NCDatasets
using Zarr
using Dates: now, format

const _OUTPUT_WRITER_TYPES = Dict{String, Type}(
    "JLD2" => JLD2Writer,
    "NetCDF" => NetCDFWriter,
    "Zarr" => ZarrWriter,
)

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
        "--free-surface-substeps", "-S"
            help = "Number of substeps to use in split-explicit free surface scheme (defaults to adaptive if not specified)"
            arg_type = Int
        "--output-directory", "-O"
            help = "Directory to write outputs to"
            arg_type = String
            default = "."
        "--cpu"
            help = "Run on CPU (rather than GPU, the default)"
            action = :store_true
        "--mpi"
            help = "Use MPI to distribute computations"
            action = :store_true
        "--ranks-along-x"
            help = "Number of ranks to distribute x dimension of grid along if using MPI"
            arg_type = Int
            default = 1
        "--use-eddy-closure"
            help = "Use a dynamic Smagorinsky eddy closure"
            action = :store_true
        "--pickup-checkpoint", "-P"
            help = "Path to checkpoint to restore simulation state from at initialisation"
            arg_type = String
        "--output-format", "-F"
            help = "Format to use for writing model outputs - one of JLD2, NetCDF or Zarr"
            arg_type = String
            default = "JLD2"
            range_tester = x -> x ∈ keys(_OUTPUT_WRITER_TYPES)
    end

    return parse_args(s)
end

function main()

    args = parse_commandline()

    args["mpi"] && MPI.Init()

    @onrank 0 @info "Parsed args\n  " * join(("$arg = $val" for (arg, val) in args), "\n  ")

    architecture = args["cpu"] ? CPU() : GPU()

    if args["mpi"]
        partition = Partition(x=args["ranks-along-x"], y=Equal())
        @onrank 0 @info partition
        architecture = Distributed(architecture; partition)
    end

    output_interval = args["output-interval-days"] * 1day

    parameters = Spall2011Parameters(;
        grid_size=Tuple(args["grid-size"]),
        surface_temperature_restoring_strength=args["surface-temperature-restoring-strength"],
        northern_basin_surface_evaporation=args["northern-basin-surface-evaporation"],
        sill_height=args["sill-height"],
        use_eddy_closure=args["use-eddy-closure"],
        split_explicit_free_surface_substeps=args["free-surface-substeps"],
    )

    mask = GyreInABox.northern_basin_mask(parameters)

    timestamp = format(now(), "yyyy-mm-dd_HH-MM-SS")

    run_directory = joinpath(args["output-directory"], "$timestamp-run")

    @onrank 0 begin
        mkdir(run_directory)
        @info "Writing outputs to $run_directory"
    end

    configuration = SimulationConfiguration(
        architecture=architecture,
        simulation_time=args["simulation-years"]*365day,
        initial_timestep=10minute,
        maximum_timestep=60minute,
        wizard_max_change=1.1,
        wizard_update_interval=10,
        output_filename=joinpath(run_directory, "spall_2012_gyre_model"),
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
        progress_message_interval=1000,
        pickup_checkpoint=args["pickup-checkpoint"],
        output_writer_type=_OUTPUT_WRITER_TYPES[args["output-format"]],
    )

    @onrank 0 begin
        jldsave(joinpath(run_directory, "parameters.jld2"); parameters)
        jldsave(joinpath(run_directory, "configuration.jld2"); configuration)
    end

    run_simulation(parameters, configuration)

end

main()
