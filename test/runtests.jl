using GyreInABox
using Oceananigans
using Test
using Aqua

@testset "GyreInABox.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GyreInABox)
    end
    @testset "Outputs" begin
        latitude_longitude_grid = LatitudeLongitudeGrid(
            size=(4, 4, 4), longitude=(0, 60), latitude=(0, 60), z=(-2000, 0), topology = (Bounded, Bounded, Bounded)
        )
        rectilinear_grid = RectilinearGrid(
            size=(4, 4, 4), x=(0, 1000), y=(0, 2000), z=(-2000, 0), topology = (Bounded, Bounded, Bounded)
        )
        all_output_types = (
            HorizontalSlice,
            LongitudeDepthSlice,
            LatitudeDepthSlice,
            XDepthSlice,
            YDepthSlice,
            DepthAveraged,
            FreeSurfaceFields,
            MOCStreamFunction,
            BarotropicStreamFunction,
        )
        grid_dependent_output_types = [
            latitude_longitude_grid => (
                HorizontalSlice,
                LongitudeDepthSlice,
                LatitudeDepthSlice,
                DepthAveraged,
                FreeSurfaceFields,
                MOCStreamFunction,
                BarotropicStreamFunction,
            ),
            rectilinear_grid => (
                HorizontalSlice,
                XDepthSlice,
                YDepthSlice,
                DepthAveraged,
                FreeSurfaceFields,
                MOCStreamFunction,
                BarotropicStreamFunction,
            ),
        ]
        @testset "Grid independent functions" begin
            for output_type in all_output_types
                output = output_type()
                @test GyreInABox.label(output) isa Symbol
                @test GyreInABox.schedule(output) isa Oceananigans.Utils.AbstractSchedule
            end
        end
        @testset "Grid dependent functions" begin
            for (grid, output_types) in grid_dependent_output_types
                for output_type in output_types
                    output = output_type()
                    model = HydrostaticFreeSurfaceModel(grid)
                    @test GyreInABox.outputs(output, model) isa NamedTuple
                    @test GyreInABox.indices(output, grid) isa Tuple
                    @test GyreInABox.axis_xlabel(output, grid) isa String
                    @test GyreInABox.axis_ylabel(output, grid) isa String
                    @test GyreInABox.axis_aspect_ratio(output, grid) isa Real
                    @test GyreInABox.axis_limits(output, grid) isa Tuple
                end
            end
        end
    end
end
