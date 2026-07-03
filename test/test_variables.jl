@testset "Variables" begin

    all_variables = (
        EastwardVelocity(),
        NorthwardVelocity(),
        UpwardVelocity(),
        Salinity(),
        Temperature(),
        FreeSurfaceDisplacement(),
        EastwardBarotropicVelocity(),
        NorthwardBarotropicVelocity(),
        MOCStreamFunction(),
        BarotropicStreamFunction(),
        NorthwardHeatTransport(),
        SpecificKineticEnergy(),
        SpecificTurbulentKineticEnergy(),
    )

    @testset "$(variable) methods" for variable in all_variables
        @test GyreInABox.short_name(variable) isa String
        @test GyreInABox.long_name(variable) isa String
        @test GyreInABox.standard_name(variable) isa String
        @test GyreInABox.units(variable) isa String
        @test GyreInABox.color_map(variable) isa Symbol
        @test GyreInABox.additional_attributes(variable) isa Dict
        @test GyreInABox.spatial_dimensions(variable) isa GyreInABox.SpatialDimensions
    end

    grid = make_test_grids()[1]

    model = make_test_model(grid)

    @testset "$(variable) field" for variable in all_variables
        @test GyreInABox.field(variable, model) isa Field
    end

end
