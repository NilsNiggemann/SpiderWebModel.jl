using Test
import SpiderWebModel as SW

@testset "StencilSpinConfig" begin
    
    S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
    oldVals = SW.getPlaquetteStencil(S,5,2)
    SW.applyPlaquette!(S,5,2,1)
    newVals = SW.getPlaquetteStencil(S,5,2)
    @test oldVals == -newVals

    S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
    oldVals = SW.getPlaquetteStencil(S,2,1)
    SW.applyPlaquette!(S,2,1,-1)
    newVals = SW.getPlaquetteStencil(S,2,1)
    @test oldVals == -newVals

end
##
Conf = SW.periodicState5x5(10)
@testset "applicablePlaquettes" begin
    @test SW.getApplicablePlaquettes(Conf) == [
        (2, 5)
        (3, 2)
        (4, 9)
        (5, 6)
        (6, 3)
        (8, 7)
        (9, 4)
    ]
end
##
H = SW.generateHilbertSpace(Conf)

eig = SW.SolveHKrylov(H.H)
##
@testset "Hamiltonian" begin # seven non-interacting plaquettes
    StatesAnalytical = Set(collect(SW.SBitVector(UInt(i)) for i = 0:(2^7-1)))

    @test Set(H.AllStates) == StatesAnalytical
    @test eig.values[1] ≈ -7
end
