using ArgParse
using GyreInABox
using JLD2
using Oceananigans
using Oceananigans.DistributedComputations
using Oceananigans.Units
using MPI
using NCDatasets
using Zarr
if Sys.isapple() && Sys.ARCH === :aarch64
    using Metal
elseif !isnothing(Sys.which("nvidia-smi"))
    using CUDA
elseif !isnothing(Sys.which("rocminfo"))
    using AMDGPU
end

const _OUTPUT_WRITER_TYPES =
    Dict{Symbol,Type}(:jld2 => JLD2Writer, :netcdf => NetCDFWriter, :zarr => ZarrWriter)

const _FLOATING_POINT_TYPES = Dict{Symbol,Type}(:float64 => Float64, :float32 => Float32)

function get_architecture(architecture_arg)
    if architecture_arg === :cpu
        CPU()
    elseif architecture_arg == :cuda
        @assert CUDA.functional()
        GPU(CUDA.CUDABackend(always_inline = true))
    elseif architecture_arg == :metal
        @assert Metal.functional()
        GPU(Metal.MetalBackend())
    elseif architecture_arg === :rocm
        @assert AMDGPU.functional()
        GPU(AMDGPU.ROCBackend())
    else
        error("Unrecognised architecture argument $architecture_arg")
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--surface-temperature-restoring-strength", "-R"
        help = "Surface temperature restoring strength / W m⁻² K⁻¹"
        arg_type = Float64
        default = 20.0
        "--northern-basin-surface-evaporation", "-E"
        help = "Reference surface net evaporation (salinity flux) in northern basin above sill / m s⁻¹"
        arg_type = Float64
        default = -2e-8
        "--northern-basin-surface-evaporation-delta", "-D"
        help = "Max change in surface salinity flux in northern basin during triangular ramp / m s⁻¹"
        arg_type = Float64
        default = -10e-8
        "--ramp-start-years"
        help = "Start time for triangular ramp of surface salinity flux in northern basin / years"
        arg_type = Float64
        default = 25.0
        "--ramp-end-years"
        help = "End time for triangular ramp of surface salinity flux in northern basin / years"
        arg_type = Float64
        default = 75.0
        "--sill-height", "-H"
        help = "Sill height / m"
        arg_type = Float64
        default = 1000.0
        "--simulation-years", "-Y"
        help = "Number of simulated year to run for"
        arg_type = Float64
        default = 100.0
        "--output-interval-days", "-I"
        help = "Interval at which to record outputs at in simulated days"
        arg_type = Float64
        default = 30.0
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
        "--architecture", "-A"
        help = "Architecture to run simulations with"
        arg_type = Symbol
        default = :cpu
        range_tester = x -> x ∈ (:cpu, :cuda, :rocm, :metal)
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
        arg_type = Symbol
        default = :jld2
        range_tester = x -> x ∈ keys(_OUTPUT_WRITER_TYPES)
        "--floating-point-type", "-T"
        help = "Floating point type to use for simulation"
        arg_type = Symbol
        default = :float64
        range_tester = x -> x ∈ keys(_FLOATING_POINT_TYPES)
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    args["mpi"] && MPI.Init()

    @onrank 0 @info "Parsed args\n  " * join(("$arg = $val" for (arg, val) in args), "\n  ")

    if args["architecture"] === :metal && args["floating-point-type"] !== :float32
        error("Metal GPUs only support float32 - use --floating-point-type float32")
    end

    architecture = get_architecture(args["architecture"])

    FT = _FLOATING_POINT_TYPES[args["floating-point-type"]]

    Oceananigans.defaults.FloatType = FT

    if args["mpi"]
        partition = Partition(; x = args["ranks-along-x"], y = Equal())
        @onrank 0 @info partition
        architecture = Distributed(architecture; partition)
    end

    output_interval = args["output-interval-days"] * 1day

    snapshot_schedule = TimeInterval(FT(output_interval))
    average_schedule =
        AveragedTimeInterval(FT(output_interval); window = FT(output_interval))

    parameters = Spall2011Parameters{FT,typeof(args["free-surface-substeps"])}(;
        grid_size = Tuple(args["grid-size"]),
        surface_temperature_restoring_strength = args["surface-temperature-restoring-strength"],
        northern_basin_surface_evaporation = args["northern-basin-surface-evaporation"],
        northern_basin_surface_evaporation_ramp_delta = args["northern-basin-surface-evaporation-delta"],
        northern_basin_surface_evaporation_ramp_start = args["ramp-start-years"]*365day,
        northern_basin_surface_evaporation_ramp_end = args["ramp-end-years"]*365day,
        sill_height = args["sill-height"],
        use_eddy_closure = args["use-eddy-closure"],
        split_explicit_free_surface_substeps = args["free-surface-substeps"],
    )

    mask = GyreInABox.northern_basin_mask(parameters)

    output_directory = args["output-directory"]

    @onrank 0 begin
        mkpath(output_directory)
        @info "Writing outputs to $output_directory"
    end

    configuration = SimulationConfiguration(;
        architecture = architecture,
        simulation_time = args["simulation-years"]*365day,
        initial_timestep = 10minute,
        maximum_timestep = 60minute,
        wizard_max_change = 1.1,
        wizard_update_interval = 10,
        output_directory = output_directory,
        output_filename_stem = "spall_2012_gyre_model",
        output_types = (
            horizontal_slice_output(; depth = 0.0, schedule = snapshot_schedule),
            horizontal_slice_output(;
                depth = (parameters.bottom_depth + parameters.sill_height),
                schedule = TimeInterval(output_interval),
            ),
            x_depth_slice_output(;
                y_or_latitude = parameters.sill_center_y,
                schedule = snapshot_schedule,
            ),
            y_depth_slice_output(;
                x_or_longitude = (parameters.domain_size_x / 2),
                schedule = snapshot_schedule,
            ),
            free_surface_output(; schedule = snapshot_schedule),
            stream_functions_output(; schedule = average_schedule),
            moc_strength_at_y_output(;
                y_or_latitude = parameters.sill_center_y,
                schedule = average_schedule,
            ),
            northward_heat_transport_at_y_output(;
                y_or_latitude = parameters.sill_center_y,
                schedule = average_schedule,
            ),
            spatially_averaged_output(;
                variables = (SpecificKineticEnergy(),),
                mask = mask,
                schedule = average_schedule,
            ),
            horizontally_averaged_output(; mask = mask, schedule = average_schedule),
        ),
        progress_message_interval = 1000,
        pickup_checkpoint = args["pickup-checkpoint"],
        output_writer_type = _OUTPUT_WRITER_TYPES[args["output-format"]],
    )

    @onrank 0 begin
        jldsave(joinpath(output_directory, "parameters.jld2"); parameters)
        jldsave(joinpath(output_directory, "configuration.jld2"); configuration)
    end

    run_simulation(parameters, configuration)
end

main()
