@testset "Spatial processors" begin
    processors = (
        GyreInABox.SpatialSliceProcessor(z=GyreInABox.SlicedSpatialDimension(0.0)),
        GyreInABox.SpatialSliceProcessor(
            x=GyreInABox.SlicedSpatialDimension(0.0),
            y=GyreInABox.SlicedSpatialDimension(0.0),
        ),
        GyreInABox.SpatialMaximumProcessor((3,)),
        GyreInABox.SpatialMaximumProcessor((1, 2)),
        GyreInABox.SpatialAverageProcessor((2,)),
        GyreInABox.SpatialAverageProcessor(
            (1, 3), GyreInABox.HorizontalCircularRegionMask(0.0, 0.0, 10.0)
        ),
        [
            GyreInABox.SpatialMaximumProcessor((2,)),
            GyreInABox.SpatialSliceProcessor(; x=GyreInABox.SlicedSpatialDimension(0.0)),
        ],
    )
    grid = make_test_grids()[1]
    field = Field{Center,Center,Center}(grid)
    dimensions = GyreInABox.SpatialDimensions(
        GyreInABox.FreeSpatialDimension(),
        GyreInABox.FreeSpatialDimension(),
        GyreInABox.FreeSpatialDimension(),
    )

    @testset "$(processor) method types" for processor in processors
        @test GyreInABox.label(processor) isa String
        @test GyreInABox.process(processor, field) isa Field
        @test GyreInABox.process(processor, dimensions) isa GyreInABox.SpatialDimensions
    end

    @testset "label" begin
        @test GyreInABox.label(
            GyreInABox.SpatialSliceProcessor(z=GyreInABox.SlicedSpatialDimension(-1.0))
        ) == "spatial_slice_at_z_-1.0"
        @test GyreInABox.label(GyreInABox.SpatialMaximumProcessor((3,))) ==
            "spatial_maximum_over_z"
        @test GyreInABox.label(GyreInABox.SpatialAverageProcessor((1, 2))) ==
            "spatial_average_over_xy"
    end

    @testset "process dimension" begin
        @test GyreInABox.process(
            GyreInABox.SpatialAverageProcessor((3,)), dimensions
        ) == GyreInABox.SpatialDimensions(
            GyreInABox.FreeSpatialDimension(),
            GyreInABox.FreeSpatialDimension(),
            GyreInABox.ReducedSpatialDimension(),
        )
        @test GyreInABox.process(
            GyreInABox.SpatialSliceProcessor(
                x=GyreInABox.SlicedSpatialDimension(0.0),
                y=GyreInABox.SlicedSpatialDimension(0.0),
            ),
            dimensions,
        ) == GyreInABox.SpatialDimensions(
            GyreInABox.SlicedSpatialDimension(0.0),
            GyreInABox.SlicedSpatialDimension(0.0),
            GyreInABox.FreeSpatialDimension(),
        )
    end
end
