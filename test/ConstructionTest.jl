using Test
using Random
import SpiderWebModel as SW
##
@testset "plaquetteList" begin
    @test length(SW.ALLGS_S12) == 70
    @test length(SW.ALLGS_S12_NOMISSING) == 140
end

Random.seed!(1234)
Conf = SW.constructConfigPath(8, 8, SW.ALLGS_S12)
Conf2 = SW.constructConfigPath(8, 8, SW.ALLGS_S12_NOMISSING)
ConfSpiral = SW.constructConfigPath(8, 8, SW.ALLGS_S12, SW.spiralPath)

@testset "ground state construction" begin
    @test SW.fulFillsConstraint(Conf)
    @test SW.fulFillsConstraint(Conf2)
    @test SW.fulFillsConstraint(ConfSpiral)
end
##
@testset "JuMP construction" begin

    L = 20
    # sols = SW.constructGroundstates(L,100,1/8)
    sols = SW.floatSpinConfig.(SW.constructGroundstates(L, 30, 1 / 7), 1 / 2)
    @test length(sols) > 0
    @test all(SW.fulFillsConstraint, sols)
    # @test all(SW.fulFillsConstraint ∘ SW.floatSpinConfig,sols)
end
