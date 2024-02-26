using SpiderWebModel, Test
import SpiderWebModel as SW
##
@testset "plaquettes" begin
    c = [1, 1, 1, 1, 1, 1, 1, 1]
    @test SW.constraint(c) == 0

    c = [1, 1, 1, -1, 1, 1, 1, 1] .* 0.5
    @test SW.constraint(c) == 1
end

##
@testset "appending path" begin
    path = Set([134, 13, 67, 87])
    SW.updatePath!(path, 341)
    @test path == Set([134, 341, 13, 67, 87])

    SW.updatePath!(path, 341)
    @test path == Set([134, 13, 67, 87])
end
include("ConstructionTest.jl")
include("FluctuationTest.jl")
