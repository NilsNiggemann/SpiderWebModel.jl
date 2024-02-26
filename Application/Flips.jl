import SpiderWebModel as SW
using HDF5
using CairoMakie
##
Confs = SW.SpinConfig.(eachslice(h5read("ConfsRaw/Confs35.h5", "Confs"), dims = 3), 1 / 2)
##
C = copy(Confs[15])
SW.plotSpinConfig(C)
SW.flipSpinsAlongDiagonal!(C, 29, 1)
SW.plotFractons(C)
