import Pkg
# Pkg.activate("../Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
meanstd(x) = (mean(x),std(x))
##

S = SW.stencilConfig(zeros(10,10),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())

μ = 0.99
CT = SW.ContinuousTimeMethod(0.2,1,(1-μ)* 0.266*prod(size(S)),SW.Hxx_RK(μ))
ψG = SW.fullVariationalFunction(S,0.1)
##
stochReconfRes = SW.stochastic_reconfiguration(S,CT,i->round(Int,8000+ 50i),ψG,10,i -> min(1,0.05 +0.2i),SW.IterativeSRSolver();Nwalkers = 6,reconfigure=false,rel_tolerance=1e-20,equilibration_steps=100,pre_equilibration_steps=10_000)
##
ψG = SW.fullVariationalFunction(S,0.1)
ψG = SW.RKFunction()
CT = SW.ContinuousTimeMethod(0.2,1,stochReconfRes.E0[end],SW.Hxx_RK(μ))
##
@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,60,2000,ψG,equilibration_steps=1000,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:12])
##

# ψG = typeof(ψG)(stochReconfRes.params)
# ψG = SW.fullVariationalFunction(S,0.1)


@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,6,8000,ψG,equilibration_steps=1000,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:12])
##
# plotEnergies(results,CT,p=100;normalize=true,dense=true,τ = 30,Emax = stochReconfRes.E0[end]+stochReconfRes.ΔE[end]/1.5,Emin = stochReconfRes.E0[end]-stochReconfRes.ΔE[end]/1.5)
plotEnergies(results,CT,p=100;normalize=true,dense=true,τ = 30,color = :red)
plotEnergies!(resultsOld,CT,p=100;normalize=true,dense=true,τ = 30)
current_figure()

##
mProj = round(Int,20 ÷ CT.τ)
BOp = SW.PlaquetteFlipOperator(S)
@time resB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,mProj,BOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])
##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj) for res in results]

BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resB,results) ])

BValsMean = mean(BVals,dims=2)[:]
BValsStd = std(BVals,dims=2)[:]

scatterlines(eachindex(BValsMean).* CT.τ,BValsMean,label="plaquette")

band!(eachindex(BValsMean).* CT.τ,BValsMean - BValsStd , BValsMean + BValsStd,color = (:black,0.2))
current_figure()
##

RBOp = SW.RandomPlaquetteFlipOperator(S)
@time resRB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,mProj,RBOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])
##
BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resRB,results) ]) ./ length(collect(SW.plaquetteIterator(S)))

BValsMean = mean(BVals,dims=2)[:]
BValsStd = std(BVals,dims=2)[:]

scatterlines!(eachindex(BValsMean).* CT.τ,BValsMean,label="plaquette",color = (:yellow))

band!(eachindex(BValsMean).* CT.τ,BValsMean - BValsStd , BValsMean + BValsStd,color = (:yellow,0.2))
current_figure()
