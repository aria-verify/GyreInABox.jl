@testset "Utils" begin
    @testset "sigmoidal_depth_profile" begin
        f_bottom = 10.0
        f_top = 20.0
        ℓ = 1.0
        z_0 = 10.0
        @test GyreInABox.sigmoidal_depth_profile(Inf, z_0, ℓ, f_bottom, f_top) ≈ f_top
        @test GyreInABox.sigmoidal_depth_profile(z_0, z_0, ℓ, f_bottom, f_top) ≈
            (f_bottom + f_top) / 2
        @test GyreInABox.sigmoidal_depth_profile(-Inf, z_0, ℓ, f_bottom, f_top) ≈ f_bottom
    end

    @testset "hyperbolically_spaced_faces" begin
        size = 10
        lower = -100.0
        upper = 0.0
        stretching_factor = 1.0
        @test GyreInABox.hyperbolically_spaced_faces(
            1, size, lower, upper, stretching_factor
        ) ≈ lower
        @test GyreInABox.hyperbolically_spaced_faces(
            size + 1, size, lower, upper, stretching_factor
        ) ≈ upper
    end

    @testset "smooth_step" begin
        @test GyreInABox.smooth_step(-1.0) == 0.0
        @test GyreInABox.smooth_step(0.0) == 0.0
        @test GyreInABox.smooth_step(0.5) == 0.5
        @test GyreInABox.smooth_step(1.0) == 1.0
        @test GyreInABox.smooth_step(2.0) == 1.0
    end

    @testset "HorizontalCircularRegionMask" begin
        x, y, r = 0.0, 0.0, 500.0
        mask = GyreInABox.HorizontalCircularRegionMask(x, y, r)
        @test summary(mask) == "circular_region_at_x_$(x)_y_$(y)_radius_$(r)"
        grid = make_test_grids()[1]
        field = Field{Center,Center,Center}(grid)
        @test mask(1, 1, 1, grid, field) == true
        @test mask(grid.Nx, grid.Ny, 1, grid, field) == false
    end

    @testset "camel_to_snake_case" begin
        @test GyreInABox.camel_to_snake_case("HelloWorld") == "hello_world"
        @test GyreInABox.camel_to_snake_case("Foo") == "foo"
        @test GyreInABox.camel_to_snake_case("HelloWorld_foo_bar") == "hello_world_foo_bar"
    end
end
