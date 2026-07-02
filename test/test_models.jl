@testset "Models" begin
    parameters = (Spall2011Parameters(), DoubleGyreParameters())

    @testset "Model interface for $(typename(p)) defaults" for p in parameters
        grid = GyreInABox.grid(p, CPU())
        @test grid isa Oceananigans.AbstractGrid
        @test GyreInABox.boundary_conditions(p) isa NamedTuple
        @test GyreInABox.forcing(p) isa NamedTuple
        @test GyreInABox.buoyancy(p) isa Oceananigans.BuoyancyFormulations.AbstractBuoyancyFormulation
        @test GyreInABox.closure(p) isa Union{Tuple, Oceananigans.TurbulenceClosures.AbstractTurbulenceClosure}
        @test GyreInABox.coriolis(p) isa Oceananigans.Coriolis.AbstractRotation
        @test GyreInABox.tracers(p) isa Tuple
        @test GyreInABox.momentum_advection(p) isa Oceananigans.Advection.AbstractAdvectionScheme
        @test GyreInABox.tracer_advection(p) isa Oceananigans.Advection.AbstractAdvectionScheme
        @test GyreInABox.free_surface(p, grid) isa Oceananigans.Models.HydrostaticFreeSurfaceModels.AbstractFreeSurface
    end

    @testset "Model set up and initialize for $(typename(p)) defaults" for p in parameters
        model = setup_model(p)
        @test model isa Oceananigans.AbstractModel
        initialize!(model, p)
    end
end
