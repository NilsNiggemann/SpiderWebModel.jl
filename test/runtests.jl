using SpiderWebModel, Test
import SpiderWebModel as SW
##
@testset "plaquettes" begin
    c = [1, 1, 1, 1, 1, 1, 1, 1]
    @test SW.constraint(c) == 0

    c = [1, 1, 1, -1, 1, 1, 1, 1] .* 0.5
    @test SW.constraint(c) == 1
end


include("ConstructionTest.jl")
include("FluctuationTest.jl")
include("StencilConfigsTest.jl")
include("GFMCTest.jl")
include("VariationalWFTest.jl")
include("IO.jl")
