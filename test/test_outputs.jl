@testset "Outputs" begin
    grids = make_test_grids()
    times = 0.0:0.1:1.0
    output_types = (
        horizontal_slice_output(; depth=0.0),
        x_depth_slice_output(; y_or_latitude=0.0),
        y_depth_slice_output(; x_or_longitude=0.0),
        depth_averaged_output(),
        free_surface_output(),
        stream_functions_output(),
        moc_strength_at_y_output(; y_or_latitude=0.0),
        northward_heat_transport_at_y_output(; y_or_latitude=0.0),
        spatially_averaged_output(variables=(GyreInABox.SpecificKineticEnergy(),)),
        horizontally_averaged_output(),
        PointModelOutput((Salinity(), Temperature()), TimeInterval(1day), [0. 0. 0.; 0. 0. -1000.])
    )
    @testset "Grid independent functions with $(output)" for output in output_types
        @test GyreInABox.label(output) isa String
        @test GyreInABox.schedule(output) isa Oceananigans.Utils.AbstractSchedule
        @test GyreInABox.variables(output) isa Tuple{Vararg{GyreInABox.AbstractModelVariable}}
        @test GyreInABox.output_filename("test", output) isa String
        @test GyreInABox.output_filename("test", output, "jld2") isa String
    end
    @testset "Grid dependent functions with $(GyreInABox.label(output)) and $(typename(grid))" for grid in grids, output in output_types
        model = make_test_model(grid)
        @test GyreInABox.outputs(output, model) isa NamedTuple
        @test GyreInABox.indices(output, grid) isa Tuple
    end
end