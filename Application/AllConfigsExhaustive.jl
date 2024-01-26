import SpiderWebModel as SW
using CairoMakie
##
using HDF5
Lx = 6
Ly = 7
@time a = SW.constructAllConfigs(Lx,Ly,SW.ALLGS_S12)
SW.plotApplPlaquettes(SW.reconstructTiling_xDirec(Lx,Ly,rand(a),SW.ALLGS_S12))
##
# aVec = stack(a,dims = 3)
# h5write("/storage/niggeni/AllConfigs.h5","AllConfigs",aVec)