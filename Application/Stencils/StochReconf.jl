import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using MKL
##
S = SW.stencilConfig(zeros(16,16),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
ψG = SW.localPlaquetteGuidingFunction(S,0.14)
# SW.get_beta_ij(ψG) .= 0.0001
nThermal = 300
##
SW.Random.seed!(1234)
stochReconfRes = SW.stochastic_reconfiguration(S,5200,ψG,3,7e-2;Nwalkers = 120,nbra = 5,rel_tolerance=1e-3,equilibration_steps=nThermal)
##
function plotVarEn(stochReconfRes)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")

    x = eachindex(stochReconfRes.E0_i)
    errorbars!(ax,x,stochReconfRes.E0_i,stochReconfRes.ΔE_i,linestyle = :solid,whiskerwidth=5)
    lines!(ax,x,stochReconfRes.E0_i)

    # lines!(ax2,x,stochReconfRes.α_i)
    fig
end
plotVarEn(stochReconfRes)
##
ψGold = SW.PlaquetteNumberGuidingFunction(0.15)
##
# using HDF5
ψGnew = SW.LocalPlaquetteGuidingFunction(stochReconfRes.params)
# ψGnew = SW.FullVariationalGuidingFunction(h5read("../../stochReconfParams.h5","L=15/params"))
##
stochReconfRes_2 = SW.stochastic_reconfiguration(S,8200,ψGnew,3,1e-4;Nwalkers = 120,nbra = 8,rel_tolerance=1e-1,equilibration_steps=nThermal,Λ=1)


ψGnew = SW.LocalPlaquetteGuidingFunction(stochReconfRes_2.params)
##
plotVarEn(stochReconfRes_2)

##
SW.Random.seed!(12322)
nBra = 5
nThermal = 3_000 ÷ nBra

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,30,500,nBra,ψGnew,5;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:32])
##
nBraOld = 5
SW.Random.seed!(1232)
# @time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,2*80,750,nBraOld,ψGold,1;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:32])
@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,30,500,nBra,ψGnew,0;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:32])
##
# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsOld,nBraOld,p=250)
plotEnergies!(results,nBra;p=250,color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
## 