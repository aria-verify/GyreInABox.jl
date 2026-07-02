"""
$(TYPEDEF)

Configuration for ocean gyre model simulation
    
$(TYPEDSIGNATURES)
    
## Details

Variables defining overall configuration of simulation such as architecture to run on,
temporal discretization and outputs to record.

$(TYPEDFIELDS)
"""
@kwdef struct SimulationConfiguration{T,A,P,W}
    "Computational architecture to run simulation on"
    architecture::A = Oceananigans.CPU()
    "Time to simulate for / s"
    simulation_time::T = 60day
    "Initial time step to use / s"
    initial_timestep::T = 10minute
    "Maximum time step to use in time step adaptation / s"
    maximum_timestep::T = 30minute
    "Directory to write output files to"
    output_directory::String = "."
    "Stem of output file names"
    output_filename_stem::String = "gyre_model"
    "Iteration interval between progress messages"
    progress_message_interval::Int = 40
    "Model variables to show statistics of in progress messages"
    progress_message_variables::Vector{<:AbstractModelVariable} = VELOCITY_AND_TRACER_VARIABLES
    "Target (advective) CFL number for time stepping wizard"
    target_cfl::T = 0.2
    "Update (iteration) interval for time stepping wizard"
    wizard_update_interval::Int = 10
    "Maximum relative time step change in each wizard update"
    wizard_max_change::T = 1.5
    "Output types to record during simulation"
    output_types::Tuple = (
        horizontal_slice_output(depth=0.0), free_surface_output(), stream_functions_output()
    )
    "Whether to checkpoint simulation at end"
    checkpoint_at_end::Bool = true
    "Optional path to checkpoint file to pickup initial simulation state from"
    pickup_checkpoint::P = nothing
    "Output writer type"
    output_writer_type::W = JLD2Writer
end

"""
Add callback for progress updates as configured in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_progress_message_callback!(
    simulation::Oceananigans.Simulation, configuration::SimulationConfiguration
)
    variables = configuration.progress_message_variables
    iteration_format_string = "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n  "
    variables_format_string = join(
        ("max(|$(short_name(v))|) = %.2e $(units(v))" for v in variables), ", "
    )
    message_string_format = Printf.Format(iteration_format_string * variables_format_string)
    progress_message(sim) = @info(
        Printf.format(
            message_string_format,
            iteration(sim),
            prettytime(sim),
            prettytime(sim.Δt),
            prettytime(sim.run_wall_time),
            (maximum(abs, field(v, sim.model)) for v in variables)...,
        ),
    )
    add_callback!(
        simulation,
        progress_message,
        IterationInterval(configuration.progress_message_interval),
    )
    nothing
end

"""
Register output writers for output types in `output_types` in `simulation` using
output filenames with directory `output_directory` and stem `output_filename_stem`
and using output format associated `output_writer_type`.

$(SIGNATURES)
"""
function add_output_writers!(
    simulation::Oceananigans.Simulation,
    output_types::Tuple,
    output_directory::String,
    output_filename_stem::String,
    output_writer_type::Type{<:Oceananigans.AbstractOutputWriter},
)
    model = simulation.model
    # Create a grid on CPU if not already to avoid issues with computing output field
    # indices from grid using scalar operations on GPU grids
    grid = (
        isa(model.grid.architecture, CPU) ? model.grid : on_architecture(CPU(), model.grid)
    )
    for output in output_types
        kwargs = if output_writer_type == NetCDFWriter
            (; output_attributes=output_attributes(variables))
        else
            (;)
        end
        simulation.output_writers[Symbol(label(output))] = output_writer_type(
            model,
            fields(output, model);
            filename=output_filename(
                output_filename_stem, output, extension(output_writer_type)
            ),
            dir=output_directory,
            indices=indices(output, grid),
            schedule=schedule(output),
            overwrite_existing=true,
            with_halos=false,
            kwargs...,
        )
    end

    nothing
end

"""
Set up simulation for model `model` with configuration `configuration`.
    
$(SIGNATURES)
"""
function setup_simulation(
    model::Oceananigans.AbstractModel, configuration::SimulationConfiguration
)
    simulation = Simulation(
        model; Δt=configuration.initial_timestep, stop_time=configuration.simulation_time
    )
    wizard = TimeStepWizard(;
        cfl=configuration.target_cfl,
        max_change=configuration.wizard_max_change,
        max_Δt=configuration.maximum_timestep,
    )
    simulation.callbacks[:wizard] = Callback(
        wizard, IterationInterval(configuration.wizard_update_interval)
    )
    add_progress_message_callback!(simulation, configuration)
    add_output_writers!(
        simulation,
        configuration.output_types,
        configuration.output_directory,
        configuration.output_filename_stem,
        configuration.output_writer_type,
    )
    simulation
end

"""
Setup and initialize model then setup and run simulation with parameters `parameters`
and configuration `configuration`.

$(SIGNATURES)
"""
function run_simulation(
    parameters::AbstractParameters, configuration::SimulationConfiguration
)
    model = setup_model(parameters; architecture=configuration.architecture)
    initialize!(model, parameters)
    simulation = setup_simulation(model, configuration)
    pickup = !isnothing(configuration.pickup_checkpoint) && configuration.pickup_checkpoint
    run!(simulation; pickup)
    if configuration.checkpoint_at_end
        filepath = joinpath(
            configuration.output_directory,
            output_filename(configuration.output_filename_stem, "checkpoint", "jld2")
        )
        Oceananigans.checkpoint(simulation; filepath)
    end
end

"""
Plot outputs of fields recorded as simulation output.

$(SIGNATURES)

## Details

Generates and record to files plots of type `plot_output_type` for model 
outputs on grid `grid` and  for a simulation configuration `configuration`. 
Keyword arguments `kwargs`  are passed through to [`plot_output()`](@ref)
and can be used to customize plot output.
"""
function plot_outputs(
    plot_output_type::AbstractPlotOutput,
    grid::AbstractGrid,
    configuration::SimulationConfiguration;
    kwargs...,
)
    for output_type in configuration.output_types
        if is_compatible(plot_output_type, output_type)
            plot_output(
                plot_output_type,
                configuration.output_directory,
                configuration.output_filename_stem,
                output_type,
                grid;
                kwargs...,
            )
        end
    end
end

function plot_outputs(
    plot_output_type::AbstractPlotOutput,
    parameters::AbstractParameters,
    configuration::SimulationConfiguration;
    kwargs...,
)
    plot_outputs(plot_output_type, grid(parameters, CPU()), configuration; kwargs...)
end

