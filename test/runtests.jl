using GyreInABox
using Oceananigans
using CairoMakie
using Test
using Aqua

@testset "GyreInABox.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GyreInABox)
    end
    @testset "Outputs" begin
        latitude_longitude_grid = LatitudeLongitudeGrid(
            size=(4, 4, 4),
            longitude=(0, 60),
            latitude=(0, 60),
            z=(-2000, 0),
            topology=(Bounded, Bounded, Bounded),
        )
        rectilinear_grid = RectilinearGrid(
            size=(4, 4, 4),
            x=(0, 1000),
            y=(0, 2000),
            z=(-2000, 0),
            topology=(Bounded, Bounded, Bounded),
        )
        times = 0.:0.1:1.
        all_output_types = (
            HorizontalSlice(),
            LongitudeDepthSlice(),
            LatitudeDepthSlice(),
            XDepthSlice(),
            YDepthSlice(),
            DepthAveraged(),
            FreeSurfaceFields(),
            MOCStreamFunction(),
            BarotropicStreamFunction(),
            MOCStrength(; y_or_latitude=1000.),
            MeridionalHeatTransport(; y_or_latitude=1000.),
            AverageKineticEnergy(),
            HorizontallyAveragedTracers(),
        )
        grid_dependent_output_types = [
            latitude_longitude_grid =>
                filter(o -> !isa(o, Union{XDepthSlice,YDepthSlice}), all_output_types),
            rectilinear_grid => filter(
                o -> !isa(o, Union{LongitudeDepthSlice,LatitudeDepthSlice}),
                all_output_types,
            ),
        ]
        @testset "Grid independent functions" begin
            for output in all_output_types
                @test GyreInABox.label(output) isa Symbol
                @test GyreInABox.schedule(output) isa Oceananigans.Utils.AbstractSchedule
            end
        end
        @testset "Grid dependent functions" begin
            for (grid, output_types) in grid_dependent_output_types
                for output in output_types
                    model = HydrostaticFreeSurfaceModel(grid; tracers=(:T, :S))
                    @test GyreInABox.outputs(output, model) isa NamedTuple
                    @test GyreInABox.indices(output, grid) isa Tuple
                    @test GyreInABox.axis_xlabel(output, grid) isa String
                    @test GyreInABox.axis_ylabel(output, grid) isa String
                    @test GyreInABox.axis_aspect_ratio(output, grid) isa Union{AxisAspect, Nothing}
                    @test GyreInABox.axis_xlimits(output, grid, times) isa Union{Tuple, Nothing}
                    @test GyreInABox.axis_ylimits(output, grid, times) isa Union{Tuple, Nothing}
                    @test GyreInABox.axis_limits(output, grid, times) isa Tuple
                end
            end
        end
    end
end
