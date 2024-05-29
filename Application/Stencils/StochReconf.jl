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
ψG = SW.constructVariationalFunction(S,0.14)
# SW.get_beta_ij(ψG) .= 0.0001
nThermal = 1000
##
res = SW.startManyWalkerGFMC(S,48,100,12,ψG,1;equilibration_steps=nThermal)
##
a = SW.reconf_obs(S,eachslice(res.SaveConfigs,dims=(3,4)),ψG,1)
##
SW.stochastic_reconfiguration(S,100,ψG,5,3e-2;Nwalkers = 10,nbra = 60,rel_tolerance=1e-1,equilibration_steps=nThermal)
##
SW.Random.seed!(1234)
stochReconfRes = SW.stochastic_reconfiguration(S,1000,ψG,20,1e-3;Nwalkers = 50,nbra = 25,rel_tolerance=1e-3,equilibration_steps=nThermal)
##
let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")

    x = eachindex(stochReconfRes.E0_i)
    errorbars!(ax,x,stochReconfRes.E0_i,stochReconfRes.ΔE_i,linestyle = :solid,whiskerwidth=5)
    lines!(ax,x,stochReconfRes.E0_i)

    # lines!(ax2,x,stochReconfRes.α_i)
    fig
end
##
ψGold = SW.PlaquetteNumberGuidingFunction(0.15)
##
# using HDF5
ψGnew = SW.VariationalGuidingFunction(stochReconfRes.params)
# ψGnew = SW.VariationalGuidingFunction(h5read("../../stochReconfParams.h5","L=15/params"))
##
stochReconfRes_2 = SW.stochastic_reconfiguration(S,3000,ψGnew,20,1e-3;Nwalkers = 50,nbra = 25,rel_tolerance=1e-3,equilibration_steps=nThermal)
##

SW.Random.seed!(1232)
nThermal = 5000

nBra = 10
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,30,120_000÷nBra,nBra,ψGnew,1;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:12])
##
nBraOld = 10
SW.Random.seed!(1232)
@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,30,120_000÷nBraOld,nBraOld,ψGold,1;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:12])

##
# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsOld,nBraOld)
plotEnergies!(results,nBra;color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
##