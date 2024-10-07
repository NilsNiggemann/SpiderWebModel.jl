import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
# using MKL
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(16,16),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
S .= SW.periodicStateDenseLoops(size(S,1))
# S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)

# ψG = SW.fullVariationalFunction(S,0.12)
ψGold = SW.PlaquetteNumberGuidingFunction(0.12)
ψG = SW.orderGuidingFunction(S,0.12)
# ψG = SW.RBM(S,1)
# SW.get_b_j(ψG) .= 0.01
# SW.get_alpha_i(ψG) .= 1
# SW.get_W_ij(ψG) .= 0.01
# SW.get_beta_ij(ψG) .= 0.0001
nThermal = 300
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.1))

##
stochReconfRes = SW.stochastic_reconfiguration(S,CT,400,ψG,100,0.3,SW.IterativeSRSolver();Nwalkers = 30,reconfigure=true,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 5,
reset = true,
# outfile = "tempSR/SR2.h5"
)
plotVarEn(stochReconfRes)
##
SW.Random.seed!(1234)

@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,30,10_000,ψGold;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:12])
##
ψGnew = deepcopy(ψG)
ψGnew.params .= stochReconfRes.params

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,30,10_000,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:12])
##

# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsOld,CT,nThermal=100,p=30,normalize=false,dense=true,τ = 10)
plotEnergies!(results,CT;nThermal=100,p=30,color=:red,normalize=false,dense=true,τ = 10) # L=15
# plotEnergies!(resultsPlaq,CT;nThermal=100,p=30,color=:blue,normalize=false,dense=true,τ = 20) # L=15
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