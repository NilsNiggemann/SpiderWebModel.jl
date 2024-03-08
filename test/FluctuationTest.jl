using Test
import SpiderWebModel as SW

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
H = SW.generateHamiltonian(Conf)

eig = SW.SolveHKrylov(H.H)
##
@testset "Hamiltonian" begin # seven non-interacting plaquettes
    StatesAnalytical = Set(collect(SW.SBitVector(UInt(i)) for i = 0:(2^7-1)))

    @test Set(H.AllStates) == StatesAnalytical
    @test eig.values[1] ≈ -7
end
