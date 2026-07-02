@testset "Dimensions" begin
    dim_types = (
        GyreInABox.FreeSpatialDimension(),
        GyreInABox.ReducedSpatialDimension(),
        GyreInABox.SlicedSpatialDimension(0.0),
    )

    grids = make_test_grids()

    times = collect(0.0:0.1:10.0)

    IndexType = Union{Colon,UnitRange,Integer}

    @testset """
        SpatialDimensions{$(typename(x_dim)), $(typename(y_dim)), $(typename(z_dim))} 
        with $(typename(grid)) method types 
        """ for x_dim in dim_types, y_dim in dim_types, z_dim in dim_types, grid in grids
        dimensions = GyreInABox.SpatialDimensions(x_dim, y_dim, z_dim)
        @test GyreInABox.indices(dimensions, grid) isa Tuple{IndexType,IndexType,IndexType}
        if !(x_dim === y_dim === z_dim === GyreInABox.FreeSpatialDimension())
            # Axis related methods intentionally not defined for dimensions with 3 free coordinates
            @test GyreInABox.axis_aspect_ratio(dimensions) isa Union{Nothing,DataAspect}
            @test GyreInABox.axis_xlabel(dimensions, grid) isa String
            @test GyreInABox.axis_ylabel(dimensions, grid) isa String
            @test GyreInABox.axis_xlimits(dimensions, grid, times) isa
                Tuple{Float64,Float64}
            @test GyreInABox.axis_ylimits(dimensions, grid, times) isa
                Union{Tuple{Float64,Float64},Nothing}
        end
    end

    @testset "update method" begin
        original_dimension = GyreInABox.FreeSpatialDimension()
        new_dimension = GyreInABox.ReducedSpatialDimension()
        dimensions = GyreInABox.SpatialDimensions(
            original_dimension, original_dimension, original_dimension
        )
        updated_dimensions = GyreInABox.update(dimensions, x=new_dimension)
        @test updated_dimensions.x == new_dimension
        @test updated_dimensions.y == original_dimension
        @test updated_dimensions.z == original_dimension
        updated_dimensions = GyreInABox.update(dimensions, y=new_dimension, z=new_dimension)
        @test updated_dimensions.x == original_dimension
        @test updated_dimensions.y == new_dimension
        @test updated_dimensions.z == new_dimension
    end
end
