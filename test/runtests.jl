using GyreInABox
using Test
using Aqua

@testset "GyreInABox.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GyreInABox)
    end
    # Write your tests here.
end
