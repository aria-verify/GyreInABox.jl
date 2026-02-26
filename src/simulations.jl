"""
$(TYPEDEF)

Configuration for ocean gyre model simulation
    
$(TYPEDSIGNATURES)
    
## Details

Variables defining overall configuration of simulation such as architecture to run on,
temporal discretization and outputs to record.

$(TYPEDFIELDS)
"""
@kwdef struct SimulationConfiguration{T,A}
    "Computational architecture to run simulation on"
    architecture::A = Oceananigans.CPU()
    "Time to simulate for / s"
    simulation_time::T = 60day
    "Initial time step to use / s"
    initial_timestep::T = 10minute
    "Maximum time step to use in time step adaptation / s"
    maximum_timestep::T = 30minute
    "Stem of output file names"
    output_filename::String = "gyre_model"
    "Iteration interval between progress messages"
    progress_message_interval::Int = 40
    "Target (advective) CFL number for time stepping wizard"
    target_cfl::T = 0.2
    "Update (iteration) interval for time stepping wizard"
    wizard_update_interval::Int = 10
    "Maximum relative time step change in each wizard update"
    wizard_max_change::T = 1.5
    "Output types to record during simulation"
    output_types::Tuple = (
        HorizontalSlice(), FreeSurfaceFields(), BarotropicStreamFunction()
    )
end

"""
Add callback for progress updates as configured in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_progress_message_callback!(
    simulation::Oceananigans.Simulation, configuration::SimulationConfiguration
)
    fields = merge(simulation.model.velocities, simulation.model.tracers)
    iteration_format_string = "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n  "
    variables_format_string = join(
        ("max(|$(variable)|) = %.2e $(unit(variable))" for variable in keys(fields)), ", "
    )
    message_string_format = Printf.Format(iteration_format_string * variables_format_string)
    progress_message(sim) = println(
        Printf.format(
            message_string_format,
            iteration(sim),
            prettytime(sim),
            prettytime(sim.Δt),
            prettytime(sim.run_wall_time),
            (maximum(abs, field) for field in values(fields))...,
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
        simulation, configuration.output_types, configuration.output_filename
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
    run!(simulation)
end

"""
Record animations of fields recorded as simulation output.

$(SIGNATURES)

## Details

Generates and record to files animations of model outputs on grid `grid` and 
for a simulation configuration `configuration`. Keyword arguments `kwargs` 
are passed through to [`record_animation()`](@ref) and can be used to customize
plot output.
"""
function record_animations(
    grid::AbstractGrid, configuration::SimulationConfiguration; kwargs...
)
    for output_type in configuration.output_types
        record_animation(configuration.output_filename, output_type, grid, kwargs...)
    end
end

function record_animations(
    parameters::AbstractParameters, configuration::SimulationConfiguration; kwargs...
)
    record_animations(grid(parameters, CPU()), configuration; kwargs...)
end

