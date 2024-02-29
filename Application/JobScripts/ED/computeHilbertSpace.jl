using Test
import SpiderWebModel as SW
cd(@__DIR__)
L = 14

# Conf = SW.getPeriodic(SW.getStairCase(L))
# tc = copy(Conf)
# tc .= -0.5

# Conf = SW.SpinConfig(SW.CircularArrays.CircularArray(Array(SW.getStairCase(L))),1/2)
Conf = SW.SpinConfig(SW.periodicPaddedMatrix(Array(SW.getStairCase(L)),1),1/2)
# SW.flipPlaquette!(tc,1,1)
# SW.updateGhostCells!(Conf.Mat)
# SW.plotApplPlaquettes(SW.SpinConfig(Conf.Mat.A_pad,1/2))
SW.plotApplPlaquettes(Conf)
##

@time H = SW.generateHilbertSpace(Conf)
# dir = mkpath("../../../Data/ED/Staircase/periodicBoundaries")
# SW.h5saveHilbertSpace(joinpath(dir,"L_$L.h5"), H)

##
L = 12
ConfOpen = SW.getStairCase(L)
@time Hopen = SW.generateHilbertSpace(ConfOpen)