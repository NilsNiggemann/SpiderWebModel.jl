import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using MKL
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(14,14),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
createLineon!(S,L) = S[end÷2+1:2:end÷2+L,end÷2+1] .= 1
createLineon!(S,1)
SW.plotFractons(S)
##
SW.Random.seed!(1234)
DT = SW.DiscreteTimeMethod(0.,3,0.2*length(S))
ψG = SW.fullVariationalFunction(S,0.2)

stochReconfRes = SW.stochastic_reconfiguration(S,DT,i->round(Int,1000+ 10*i),ψG,40,i->5*max(2. - 1*log(i),0.1) ,SW.IterativeSRSolver();Nwalkers = 6*30,reconfigure=true,rel_tolerance=1e-8,equilibration_steps=100,pre_equilibration_steps=50_000,scatter_fraction=0.6)
plotVarEn(stochReconfRes)
##
DT = SW.DiscreteTimeMethod(0.,7,0.2*length(S))

ψG = SW.FullVariationalGuidingFunction(stochReconfRes.params_steps[:,:,end])
##
scatter_fraction = 0.2
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,320,2500,ψG;equilibration_steps=500,pre_equilibration_steps=50_000,scatter_fraction) for i in 1:6])
##
plotEnergies(results,DT,p=50)
##
S = SW.stencilConfig(zeros(14,14),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
createLineon!(S,L) = S[end÷2+1:2:end÷2+L,end÷2+1] .= 1
createLineon!(S,5)
SW.plotFractons(S)
##
SW.Random.seed!(1234)
DT = SW.DiscreteTimeMethod(0.,3,0.2*length(S))
ψG = SW.fullVariationalFunction(S,0.2)

stochReconfRes5 = SW.stochastic_reconfiguration(S,DT,i->round(Int,1000+ 10*i),ψG,40,i->5*max(2. - 1*log(i),0.1) ,SW.IterativeSRSolver();Nwalkers = 6*30,reconfigure=true,rel_tolerance=1e-8,equilibration_steps=100,pre_equilibration_steps=50_000,scatter_fraction=0.6)
plotVarEn(stochReconfRes5)
##
DT = SW.DiscreteTimeMethod(0.,7,0.2*length(S))

ψG = SW.FullVariationalGuidingFunction(stochReconfRes.params_steps[:,:,end])
##
scatter_fraction = 0.2
@time results5 = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,320,2000,ψG;equilibration_steps=500,pre_equilibration_steps=50_000,scatter_fraction) for i in 1:6])
##
plotEnergies(results,DT,p=120)
plotEnergies!(results5,DT,p=120,color = :red)
current_figure()