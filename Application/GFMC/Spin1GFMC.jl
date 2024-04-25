cd(@__DIR__)

L = 40
nBra = 10
NSteps = 20_000
NWalkers = 100
α = 0.15
λ = 1

outfile = "../Data/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_alpha=$(α)_Lam=$λ.h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5, H5Zblosc
##
h5open(outfile,"w") do f
    f["test",blosc=9] = rand(Int8,20,20)
end
rm(outfile)
#___________Spin-1_______________________
##
parentState = SW.stencilConfig(zeros(L,L),1)

GuidingWaveFunction(x) = SW.varitationalFunc(α,x,0)
##
@info "starting run"  
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,varFuncTest,1) for _ in 1:8])
@time results = SW.startManyWalkerGFMC(parentState,NWalkers,NSteps÷nBra,nBra,GuidingWaveFunction,λ)
##
h5open(outfile,"w") do f
    f["TotalWeights",blosc=9] = results.TotalWeights
    f["energies",blosc=9] = results.energies
    f["SaveConfigs",blosc=9] = results.SaveConfigs
    f["reconfigurationTable",blosc=9] = results.reconfTable
end