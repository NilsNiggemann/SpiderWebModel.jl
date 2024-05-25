import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

##
S = SW.stencilConfig(zeros(12,12),1;
# boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
ψG = SW.constructVariationalFunction(S,0.15)
# SW.get_beta_ij(ψG) .= 0.01
nThermal = 100

##
SW.repeatStochReconf(S,200,ψG,5,5e-2;Nwalkers = 12,nbra = 10,error_threshold=1e-1,equilibration_steps=nThermal)
##
stochReconfRes = SW.repeatStochReconf(S,200,ψG,10,5e-2;Nwalkers = 12,nbra = 10,error_threshold=1e-1,equilibration_steps=nThermal)
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
# ψGnew = SW.VariationalGuidingFunction(stochReconfRes.params)
ψGold = SW.PlaquetteNumberGuidingFunction(0.15)
##
nBra = 6
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,5,50_000÷nBra,nBra,ψGnew,1;equilibration_steps=nThermal) for _ in 1:24])

##
# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8)
plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
