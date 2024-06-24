import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using MKL
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(18,18),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
# S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)

ψG = SW.fullVariationalFunction(S,0.15)
# ψG = SW.localPlaquetteGuidingFunction(S,0.15)
# SW.get_beta_ij(ψG) .= 0.0001
nThermal = 300
##
SW.Random.seed!(1234)
DT = SW.DiscreteTimeMethod(0.,10,prod(size(S)))

stochReconfRes = SW.stochastic_reconfiguration(S,DT,i->round(Int,10+ 2*i),ψG,50,0.03,SW.IterativeSRSolver();Nwalkers = 6*20,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000)
##
function plotVarEn(stochReconfRes)
    fig = Figure(theme = theme_SimpleTicks())

    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Energy",xlabelvisible=false,xticklabelsvisible=false)
    ax2 = Axis(fig[2,1], xlabel = "Iteration", ylabel = "α")

    E_min = minimum(stochReconfRes.E0)
    Epl =stochReconfRes.E0
    # Epl = abs.(stochReconfRes.E0_i .- E_min)
    x = eachindex(Epl)
    errorbars!(ax,x,Epl,stochReconfRes.ΔE,whiskerwidth=5)
    lines!(ax,x,Epl)
    # ylims!(ax,0.001,10.)
    # lines!(ax2,x,SW.mean(stochReconfRes.params_steps,dims=1)[1,1,:])
    lines!(ax2,x,stochReconfRes.params_steps[1,1,:])
    # lines!(ax2,x,SW.mean(stochReconfRes.params_steps,dims=1)[1,end,:])
    fig
end
plotVarEn(stochReconfRes)
##
SW.Random.seed!(12322)
nThermal = 0
DT = SW.DiscreteTimeMethod(0.,5,prod(size(S))/4)
ψGnew = typeof(ψG)(stochReconfRes.params)

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,80,1_000,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:6])
##
SW.Random.seed!(1232)

# @time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,2*80,750,nBraOld,ψGold,1;equilibration_steps=nThermal,pre_equilibration_steps=nBra*nThermal,w_avg_estimate = 8.) for _ in 1:32])
@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,80,1_000,ψG;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:6])
##
# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsOld,DT.nBranch,nThermal=100,p=400,normalize=true)
plotEnergies!(results,DT.nBranch;nThermal=100,p=400,color=:red,normalize=true) # L=15
# plotEnergies(results,DT.nBranch;nThermal=1,p=1000,color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
## 

SqsGFMC = SW.getSqsGFMC(results,350,DT.nBranch)
##
function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end

with_theme(theme_PiTicks()) do 
    Sq = mean(SqsGFMC) ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig = Figure(fontsize = 22,size = (800,400))
    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    err = sqrt.(var(SqsGFMC)) ./4
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
    hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis)

    Colorbar(fig[2,1],hmMC,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,2],hmerr,label = L"\sigma( \mathcal{S}^{zz}(\textbf{q}))",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,3],hmFT,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.8),vertical=false,ticks = SimpleTicks())

    rowsize!(fig.layout,2,Relative(0.1))
    fig
end

