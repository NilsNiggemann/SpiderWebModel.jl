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
