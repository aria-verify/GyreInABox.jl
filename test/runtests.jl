using GyreInABox
using Oceananigans
using Oceananigans.Units
using CairoMakie
using Test
using Aqua

function make_test_grids(size=(4, 4, 4), topology=(Bounded, Bounded, Bounded))
    rectilinear_grid = RectilinearGrid(
        size=size,
        x=(0, 1000),
        y=(0, 2000),
        z=(-2000, 0),
        topology=topology,
    )
    latitude_longitude_grid = LatitudeLongitudeGrid(
        size=size,
        longitude=(0, 60),
        latitude=(0, 60),
        z=(-2000, 0),
        topology=topology,
    )
    immersed_boundary_grid = ImmersedBoundaryGrid(
        rectilinear_grid, GridFittedBottom((x, y) -> -2000.0 + x + 0.5 * y)
    )
    (rectilinear_grid, latitude_longitude_grid, immersed_boundary_grid)
end

function make_test_model(grid)
    # Explicitly include T, S and e tracers and use split explicit free surface
    # to ensure all variable types have corresponding fields defined
    HydrostaticFreeSurfaceModel(
        grid; tracers=(:T, :S, :e), free_surface=SplitExplicitFreeSurface(grid)
    )
end

typename(t) = nameof(typeof(t))

@testset "GyreInABox.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GyreInABox)
    end
    include("test_variables.jl")
    include("test_dimensions.jl")
    include("test_processors.jl")
    include("test_outputs.jl")
    include("test_utils.jl")
    include("test_models.jl")
    include("test_simulations.jl")
end
