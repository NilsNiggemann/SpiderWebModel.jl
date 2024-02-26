import SpiderWebModel as SW
using Test
##

Conf = SW.periodicState5x5(10)

H = SW.generateHilbertSpace(Conf)

eig = SW.SolveHKrylov(H.H)

##
@testset "write and read Hilbert space" begin
    tmp = tempname()
    SW.h5write(tmp, H.plaqMapping)
    PM = SW.h5readPlaqMapping(tmp)
    @test PM == H.plaqMapping
end
