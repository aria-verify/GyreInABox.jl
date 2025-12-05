"""
Add callback for progress updates as configured in `configuration` in `simulation`.

$(SIGNATURES)
"""
function add_progress_message_callback!(
    simulation::Oceananigans.Simulation,
    configuration::GyreInABoxConfiguration{T}
) where {T}
    fields = merge(simulation.model.velocities, simulation.model.tracers)
    iteration_format_string = "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n  "
    variables_format_string = join(
        ("max(|$(variable)|) = %.2e $(unit(variable))" for variable in keys(fields)),
        ", "
    )
    message_string_format = Printf.Format(
        iteration_format_string * variables_format_string
    )
    progress_message(sim) = println(
        Printf.format(
            message_string_format,
            iteration(sim),
            prettytime(sim),
            prettytime(sim.Δt),
            prettytime(sim.run_wall_time),
            (maximum(abs, field) for field in values(fields))...
        )
    )
    add_callback!(
        simulation,
        progress_message,
        IterationInterval(configuration.progress_message_interval)
    )
    nothing
end

"""
Set up simulation for model `model` with configuration `configuration`.
    
$(SIGNATURES)
"""
function setup_simulation(
    model::Oceananigans.AbstractModel,
    configuration::GyreInABoxConfiguration{T}
) where {T}
    simulation = Simulation(
        model,
        Δt=configuration.initial_timestep,
        stop_time=configuration.simulation_time
    )
    wizard = TimeStepWizard(
        cfl=configuration.target_cfl,
        max_change=configuration.wizard_max_change,
        max_Δt=configuration.maximum_timestep
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
    parameters::GyreInABoxParameters{T},
    configuration::GyreInABoxConfiguration{T}
) where {T}
    model = setup_model(parameters, configuration)
    initialize_model!(model, parameters)
    simulation = setup_simulation(model, configuration)
    run!(simulation)
end

"""
Record animations of fields recorded as simulation output.

$(SIGNATURES)

## Details

Generates and record to files animations of model outputs for a model configuration
`configuration`. Keyword arguments `kwargs` are passed through to 
[`record_animation()`](@ref) and can be used to customize plot output.
"""
function record_animations(configuration::GyreInABoxConfiguration; kwargs...)
    grid = setup_grid(configuration, architecture=CPU())
    for output_type in configuration.output_types
        record_animation(configuration.output_filename, output_type, grid, kwargs...)
    end
end
