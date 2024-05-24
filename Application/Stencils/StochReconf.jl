import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

##
S = SW.stencilConfig(zeros(26,26),1;
# boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
ψG = SW.PlaquetteNumberGuidingFunction(0.05)
nThermal = 500
##

stochReconfRes = SW.repeatStochReconf(S,1000,ψG,20,1e-2;Nwalkers = 6,nbra = 5,error_threshold=1e-2,equilibration_steps=nThermal)

##
let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")

    x = eachindex(stochReconfRes.E0_i)
    errorbars!(ax,x,stochReconfRes.E0_i,stochReconfRes.ΔE_i,linestyle = :solid,whiskerwidth=5)
    lines!(ax,x,stochReconfRes.E0_i)

    lines!(ax2,x,stochReconfRes.α_i)
    fig
end