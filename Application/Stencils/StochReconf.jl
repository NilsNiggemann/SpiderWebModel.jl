import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
# using MKL
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(8,8),1;
boundaryCondition = :open
)
S .= SW.periodicStateDenseLoops(size(S,1))
# S = SW.stencilConfig(SW.getStairCase(20),1/2;
# boundaryCondition = :periodic
# )
# S = SW.stencilConfig(SW.getStairCase(12),1/2;
# boundaryCondition = :open,
# )

# ψG = SW.fullVariationalFunction(S,0.04)
# ψGold = SW.PlaquetteNumberGuidingFunction(0.12)
# ψG = SW.orderGuidingFunction(S,0.12)
# ψG.M_i .= 1e-6
# ψG = SW.localPlaquetteGuidingFunction(S,0.001)
ψG = SW.RBMSpin1(S,1,Float64)
ψGSymm = SW.getNonSymmetric(ψG)
SW.Random.seed!(1234)
# ψGSymm = SW.symmetrize(S,ψG,(4,4))
# ψG = SW.RBM(S,1)
SW.rand!(SW.get_params(ψG))
SW.get_params(ψG) .*= 1e-3
# ψGSymm = SW.symmetrize(S,ψG,(4,4))
ψGold = SW.PlaquetteNumberGuidingFunction(0.12)
nThermal = 1000
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.5))
CTSR = SW.ContinuousTimeMethod(4,Hxx = CT.Hxx)

##
cappedGrowth(x,start,stop,offset,growth) = start + 0.5(stop-start)* (1 +tanh(growth *(x -offset)))
numSteps(i) = round(Int,cappedGrowth(i,10,30,100,0.02))
learningRate(i) = cappedGrowth(i,0.01,0.005,150,0.02)
SRSteps = 300
lines(1:SRSteps,learningRate.(1:SRSteps))
lines!(1:SRSteps,numSteps.(1:SRSteps))
current_figure()
##
learningRate(i) = i > 450 ? 1e-4 : 1e-4
numSteps(i) = i > 500 ? 200 : 10
##

stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,10 ,ψGSymm,SRSteps,learningRate,SW.IterativeSRSolver();Nwalkers = 1*28,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 10,
reset = false,
# outfile = "tempSR/SR2.h5"
)

ψGnew = deepcopy(ψG)
SW.get_params(ψGnew) .= stochReconfRes.params

plotVarEn(stochReconfRes)
##

SW.Random.seed!(1234)
NWalkers = 12
NSteps = 2000
@time resultsNaive = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGold;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])

@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψG;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])
##
SW.Random.seed!(1234)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])
##

# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsNaive,CT,normalize=false,dense=true,τ = 20)
# plotEnergies!(resultsOld,CT,normalize=false,dense=true,τ = 20,color = :blue)
plotEnergies!(results,CT;color=:red,normalize=false,dense=true,τ = 20) # L=15
# plotEnergies!(resultsPlaq,CT;nThermal=100,p=30,color=:blue,normalize=false,dense=true,τ = 20) # L=15
# plotEnergies(results,DT.nBranch;nThermal=1,p=1000,color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
## 
SqsGFMC = SW.getSqsGFMC(resultsOld,1:30)
##
with_theme(theme_SimpleTicks()) do 
    Sqmean = dropmean(SqsGFMC,dims = 4)
    Sqerr = dropstd(SqsGFMC,dims = 4)

    inds = [
        (4,5),
        (5,5),
        (6,5),
        (5,4),
        (5,6),
        Tuple(argmax(Sqmean[:,:,end]))
    ]
    fig = Figure(fontsize = 22,size = (800,400))
    ax = Axis(fig[1,1],xlabel = L"\tau",ylabel = L"\mathcal{S}^{zz}(\textbf{q})")
    x = axes(SqsGFMC)[3] .* CT.τ
    for (i,j) in inds
        sqm = Sqmean[i,j,:]
        sqe = Sqerr[i,j,:]
        
        l = lines!(x,sqm)
        band!(x,sqm .- sqe,sqm .+ sqe ,color = (l.color[],0.2))

        fig
    end

    current_figure()
end
##
with_theme(theme_PiTicks()) do 
    Sq = dropmean(SqsGFMC,dims=4)[:,:,end] ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig = Figure(fontsize = 22,size = (800,400))
    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    err = dropstd(SqsGFMC,dims=4)[:,:,end] ./4
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