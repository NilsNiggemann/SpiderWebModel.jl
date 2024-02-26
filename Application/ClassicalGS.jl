import SpiderWebModel as SW
##
Lx = 6
Ly = 6
@time a = SW.constructAllConfigs(Lx, Ly, SW.ALLGS_S12)
a_rec = SW.fillEmptyStates(a, Lx, Ly, SW.ALLGS_S12)
##

@profview GS = SW.constructGSFromTiles_Threads(a_rec, 9, 9, numTries = 10000000)
##
SW.plotApplPlaquettes(rand(GS))
