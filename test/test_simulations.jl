@testset "Simulations" begin
    parameters = (
        Spall2011Parameters(; grid_size=(25, 50, 10)),
        DoubleGyreParameters(; grid_size=(30, 30, 10)),
    )

    @testset "Run simulation with $(typename(p)) with small grid" for p in parameters
        mktempdir() do tmp_dir
            configuration = SimulationConfiguration(;
                simulation_time=2day,
                output_directory=tmp_dir,
                checkpoint_at_end=true,
                initial_timestep=15minute,
                maximum_timestep=120minute,
            )
            run_simulation(p, configuration)

            @test isfile(
                joinpath(
                    tmp_dir,
                    GyreInABox.output_filename(
                        configuration.output_filename_stem, "checkpoint", "jld2"
                    ),
                ),
            )
            for output_type in configuration.output_types
                @test isfile(
                    joinpath(
                        tmp_dir,
                        GyreInABox.output_filename(
                            configuration.output_filename_stem,
                            output_type,
                            GyreInABox.extension(configuration.output_writer_type),
                        ),
                    ),
                )
            end
        end
    end
end
