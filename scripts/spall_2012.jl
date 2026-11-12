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
        "--averaged-output-interval-days", "-I"
        help = "Interval at which to record time-averaged outputs at in simulated days"
        arg_type = Float64
        default = 30.0
        "--slice-output-interval-days", "-J"
        help = "Interval at which to record snapshotted slice outputs at in simulated days"
        arg_type = Float64
        default = 30.0
        "--profile-output-interval-days", "-K"
        help = "Interval at which to record snapshotted profile outputs at in simulated days"
        arg_type = Float64
        default = 1.0
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
        "--output-types"
        help = "Output type groups to include in simulation output"
        arg_type = Symbol
        nargs = '*'
        default = [
            :scalar,
            :depth_profile,
            :surface_slice,
            :central_slice,
            :free_surface,
            :stream_function,
        ]
        "--initial-timestep-minutes"
        help = "Initial simulation timestep / minutes"
        arg_type = Float64
        default = 10.0
        "--maximum-timestep-minutes"
        help = "Maximum simulation timestep / minutes"
        arg_type = Float64
        default = 60.0
        "--timestep-wizard-update-interval"
        help = "Number of iterations between each adaptive timestep update"
        arg_type = Int
        default = 10
        "--timestep-wizard-max-change"
        help = "Maximum relative increase (multiplier) for timestep in adaptation"
        arg_type = Float64
        default = 1.1
        "--progress-message-interval"
        help = "Number of iterations between each progress message update"
        arg_type = Int
        default = 1000

    end

    return parse_args(s)
end

function profile_points(parameters)
    [
        (x = 0.05 * parameters.domain_size_x, y = parameters.sill_center_y),
        (x = 0.95 * parameters.domain_size_x, y = parameters.sill_center_y),
        (
            x = parameters.domain_size_x / 2,
            y = parameters.domain_size_y - parameters.domain_size_x / 2,
        ),
        (x = parameters.domain_size_x / 2, y = 0.975 * parameters.domain_size_y),
    ]
end

function construct_output_types(
    averaged_output_interval_days::Float64,
    slice_output_interval_days::Float64,
    profile_output_interval_days::Float64,
    output_type_groups::Vector{Symbol},
    FT::DataType,
    parameters::Spall2011Parameters,
)
    slice_schedule = TimeInterval(FT(slice_output_interval_days * day))
    profile_schedule = TimeInterval(FT(profile_output_interval_days * day))
    average_schedule = AveragedTimeInterval(
        FT(averaged_output_interval_days * day);
        window = FT(averaged_output_interval_days * day),
    )
    mask = GyreInABox.northern_basin_mask(parameters)
    output_types = GyreInABox.FieldModelOutput[]
    (:scalar in output_type_groups) && append!(
        output_types,
        (
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
        ),
    )
    (:depth_profile in output_type_groups) && append!(
        output_types,
        (
            horizontally_averaged_output(; mask = mask, schedule = average_schedule),
            (
                depth_profile_output(;
                    x_or_longitude = p.x,
                    y_or_latitude = p.y,
                    schedule = profile_schedule,
                ) for p in profile_points(parameters)
            )...,
        ),
    )
    (:surface_slice in output_type_groups) && append!(
        output_types,
        (horizontal_slice_output(; depth = 0.0, schedule = slice_schedule),),
    )
    (:central_slice in output_type_groups) && append!(
        output_types,
        (
            horizontal_slice_output(;
                depth = (parameters.bottom_depth + parameters.sill_height),
                schedule = TimeInterval(slice_schedule),
            ),
            x_depth_slice_output(;
                y_or_latitude = parameters.sill_center_y,
                schedule = slice_schedule,
            ),
            y_depth_slice_output(;
                x_or_longitude = (parameters.domain_size_x / 2),
                schedule = slice_schedule,
            ),
        ),
    )
    (:free_surface in output_type_groups) &&
        append!(output_types, (free_surface_output(; schedule = slice_schedule),))
    (:stream_function in output_type_groups) &&
        append!(output_types, (stream_functions_output(; schedule = average_schedule),))
    return output_types
end

construct_output_types(args, FT::DataType, parameters::Spall2011Parameters) =
    construct_output_types(
        args["averaged-output-interval-days"],
        args["slice-output-interval-days"],
        args["profile-output-interval-days"],
        args["output-types"],
        FT,
        parameters,
    )

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

    output_directory = args["output-directory"]

    @onrank 0 begin
        mkpath(output_directory)
        @info "Writing outputs to $output_directory"
    end

    configuration = SimulationConfiguration(;
        architecture = architecture,
        simulation_time = args["simulation-years"]*365day,
        initial_timestep = args["initial-timestep-minutes"] * minute,
        maximum_timestep = args["maximum-timestep-minutes"] * minute,
        wizard_max_change = args["timestep-wizard-max-change"],
        wizard_update_interval = args["timestep-wizard-update-interval"],
        output_directory = output_directory,
        output_filename_stem = "spall_2012_gyre_model",
        output_types = construct_output_types(args, FT, parameters),
        progress_message_interval = args["progress-message-interval"],
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
