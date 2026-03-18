using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

using ThreadPinning
ThreadPinning.pinthreads(:cores)

import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

##

# res = [SW.startManyWalkerGFMC(S,CT,100,9000,ψG;initializer) for _ in 1:6];

# Sqs = [O.StructureFactor for O in ObsRuns[findall(x->minimum(x)>-90,energies)]]
# energies = [O.Energy for O in ObsRuns[findall(x->minimum(x)>-90,energies)]]

##
function connectivity(confs,Spin)
    Lx,Ly = size(confs[1])
    degrees = zeros(Int, length(confs))
    buffer = zeros(Int, Lx*Ly)
    @info "" size(degrees) size(buffer) size(confs)

    S = SW.stencilConfig(Spin*ones(Lx,Ly),Spin,boundaryCondition= :periodic)

    for i in eachindex(degrees,confs)
        S .= confs[i]
        degrees[i] = sum(SW.getNPlaq!(buffer,S))
    end
    return degrees
end

function getAllConnectivities(SaveConfs,Spin)
    Lx,Ly = size(SaveConfs)[1:2]
    Confs_reshape = reshape(eachslice(SaveConfs,dims=(3,4)),:)
    return connectivity(Confs_reshape,Spin)
end


##
Nwalkers = 28*6
Nconfs = 2000
L = 28
# S = SW.stencilConfig(zeros(12,12),1,boundaryCondition= :periodic) .= 2SW.periodicState6x6Condensate(12)

S = SW.get4x4PeriodicSpinConf(L,6)
ψG = SW.RKFunction()
# SW.get_params(ψG) .= h5read("../../Data/GWF6x6/6x6_GF_24_mu0.h5","params")
initializer = SW.UnguidedWalkInitializer(4_000_000,1.0)

CT = SW.ContinuousTimeMethod(5,w_avg_estimate =0. ,Hxx = SW.Hxx_RK(1.0))

##
SHalf = SW.stencilConfig(0.5*ones(L,L),0.5,boundaryCondition= :periodic) .= 2SW.getStairCase(L)
initializer = SW.UnguidedWalkInitializer(4_000_000,1.0)

res_task = Threads.@spawn SW.startManyWalkerGFMC(S,CT,Nwalkers,Nconfs,ψG;initializer)
resHalf_task = Threads.@spawn SW.startManyWalkerGFMC(SHalf,CT,Nwalkers,Nconfs,ψG;initializer)

res = fetch(res_task)
resHalf = fetch(resHalf_task)
println("Done with GFMC runs")
##

conns_spin1 = Int16.(getAllConnectivities(res.SaveConfigs,1))
conns_spinhalf = Int16.(getAllConnectivities(resHalf.SaveConfigs,0.5))

##
using SpiderWebModel.HDF5
h5write("Application/Data/Connectivity/sector_analysis_28x28_staircase.h5","connectivities_spin1",conns_spin1)
h5write("Application/Data/Connectivity/sector_analysis_28x28_staircase.h5","connectivities_spinhalf",conns_spinhalf)
