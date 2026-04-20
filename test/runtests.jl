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
        immersed_boundary_grid = ImmersedBoundaryGrid(
            rectilinear_grid, GridFittedBottom((x, y) -> -2000. + x + 0.5 * y)
        )
        grids = (latitude_longitude_grid, rectilinear_grid, immersed_boundary_grid)
        times = 0.:0.1:1.
        all_output_types = (
            HorizontalSlice(),
            XDepthSlice(; y_or_latitude=0.),
            YDepthSlice(; x_or_longitude=0.),
            DepthAveraged(),
            FreeSurfaceFields(),
            MOCStreamFunction(),
            BarotropicStreamFunction(),
            MOCStrength(; y_or_latitude=0.),
            MeridionalHeatTransport(; y_or_latitude=0.),
            AverageKineticEnergy(),
            HorizontallyAveragedTracers(),
        )
        @testset "Grid independent functions" begin
            for output in all_output_types
                @test GyreInABox.label(output) isa Symbol
                @test GyreInABox.schedule(output) isa Oceananigans.Utils.AbstractSchedule
            end
        end
        @testset "Grid dependent functions" begin
            for grid in grids
                for output in all_output_types
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
