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
boundaryCondition = :periodic
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
# ψG = SW.RBMSpin1(S,1,Float64)
ψG = SW.SimpleJastrowFunction(S,Float64)
# ψG = SW.orderGuidingFunction(S)
# ψG = SW.PlaquetteRBM(S,1,Float64)
# Symmetry = SW.SymmetryGroup(SW.TranslationalSymmetry(SA[2,2],SA[-2,2]),SW.ExchangeSymmetry())
Symmetry = SW.TranslationalSymmetry(SA[2,2],SA[-2,2])
# Symmetry = SW.SymmetryGroup(SW.ExchangeSymmetry())
# ψGSymm = SW.getNonSymmetric(ψG)
ψGSymm = SW.symmetrize(ψG,Symmetry,S)
# SW.reduceParams!(ψGSymm,,S)
SW.Random.seed!(1234)
# ψGSymm = SW.symmetrize(S,ψG,(4,4))
# ψG = SW.RBM(S,1)
# SW.rand!(SW.get_params(ψG)) .*= 1e-3
SW.rand!(ψGSymm,1e-3)
ψG.v_ij .= SW.Symmetric(ψG.v_ij)
# ψG.α .= 0.1
##

# ψGSymm = SW.symmetrize(S,ψG,(4,4))
ψGold = SW.PlaquetteNumberGuidingFunction(0.12)
nThermal = 400
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.4))
CTSR = SW.ContinuousTimeMethod(15,Hxx = CT.Hxx)


##
cappedGrowth(x,start,stop,offset,growth) = start + 0.5(stop-start)* (1 +tanh(growth *(x -offset)))
SRSteps = 1500

numSteps(i) = round(Int,cappedGrowth(i,20,40,SRSteps - SRSteps÷2,10/SRSteps))
learningRate(i) = cappedGrowth(i,8e-1,4e-1,SRSteps - SRSteps÷2,5/SRSteps)
lines(1:SRSteps,5000*learningRate.(1:SRSteps))
lines!(1:SRSteps,numSteps.(1:SRSteps))
# current_figure()
##
learningRate(i) = i > SRSteps ÷8 ? 1e-2 : 1e-4
numSteps(i) = i > SRSteps ÷2 ? 50 : 10
##

stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,20 ,ψG,SRSteps,5e-4,SW.IterativeSRSolver();Nwalkers = 3*20,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 20,
reset = false,
# outfile = "tempSR/SR2.h5"
)

ψGnew = deepcopy(ψG)
SW.get_params(ψGnew) .= stochReconfRes.params
# SW.get_params(ψGnew) .= stochReconfRes.params_steps[:,begin]
plotVarEn(stochReconfRes,movavg = 30,alpha_index = 1)
##

stochReconfResSymm = SW.stochastic_reconfiguration(S,CTSR,15 ,ψGSymm,SRSteps,1e-3,SW.IterativeSRSolver();Nwalkers = 1*20,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 20,
reset = false,
# outfile = "tempSR/SR2.h5"
)

psiSymm = deepcopy(ψG)
SW.get_params(psiSymm) .= stochReconfResSymm.params
plotVarEn(stochReconfResSymm,movavg = 30,alpha_index = 1)

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22,size = (1000,600))

    axSR1 = Axis(fig[1,1],xlabel = L"iter$$",ylabel = L"E", xlabelvisible = false,xticklabelsvisible = false)
    axSR2 = Axis(fig[2,1],xlabel = L"iter$$",ylabel = L"\Delta E")
    
    errlines!(axSR1,stochReconfRes.E0,stochReconfRes.ΔE;label = "Non-symmetric",color = :black)
    errlines!(axSR1,stochReconfResSymm.E0,stochReconfResSymm.ΔE;label = "Symmetric",color = :red)
    lines!(axSR2,stochReconfRes.ΔE,label = "Non-symmetric",color = :black)
    lines!(axSR2,stochReconfResSymm.ΔE,label = "Non-symmetric",color = :red)
    axmag1 = Axis(fig[1,2],xlabel = "x",ylabel = "y",title = L"$m$ Non-symmetric",aspect=1,xlabelvisible=false,xticklabelsvisible=false)
    axmag2 = Axis(fig[1,3],xlabel = "x",ylabel = "y",title = L"$m$ Symmetric",aspect=1, ylabelvisible = false,yticklabelsvisible = false,xlabelvisible=false,xticklabelsvisible=false)
    
    mag = reshape(ψGnew.m_i,size(S))
    # magSymm = reshape(ψG.m_i,size(S))
    magSymm = reshape(psiSymm.m_i,size(S))

    cmap = min(minimum(mag),minimum(magSymm)),max(maximum(mag),maximum(magSymm))

    heatmap!(axmag1,mag;colormap = :viridis,colorrange = cmap)
    hm = heatmap!(axmag2,magSymm;colormap = :viridis,colorrange = cmap)
    Colorbar(fig[1,4],hm)

    ax = Axis(fig[2,2],xlabel = "x",ylabel = "y",title = L"$v_{ij}$ Non-symmetric",aspect=1)
    ax2 = Axis(fig[2,3],xlabel = "x",ylabel = "y",title = L"$v_{ij}$ Symmetric",aspect=1, ylabelvisible = false,yticklabelsvisible = false)

    vijFull = reshape(ψGnew.v_ij,size(S,1),size(S,2),size(S,1),size(S,2))

    # x,y = 4,3
    x,y = 4,4

    vijSymm = reshape(psiSymm.v_ij,size(S,1),size(S,2),size(S,1),size(S,2))

    v = vijFull[x,y,:,:]
    vSymm = vijSymm[x,y,:,:]

    cmap = min(minimum(v),minimum(vSymm)),max(maximum(v),maximum(vSymm))

    heatmap!(ax,v;colormap = :viridis,colorrange = cmap)
    hm = heatmap!(ax2,vSymm;colormap = :viridis,colorrange = cmap)
    scatter!(ax,[x],[y],color = :red)
    scatter!(ax2,[x],[y],color = :red)

    Colorbar(fig[2,4],hm)

    colsize!(fig.layout,1,Relative(0.4))
    axislegend(axSR1)

    fig
end

##

SW.Random.seed!(1234)
NWalkers = 20*3
NSteps = 2000
@time resultsNaive = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGold;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:20*4])

# @time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψG;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:20*4])
##
@time resultsHighAcc = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20NWalkers,2NSteps,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:10])

##
SW.Random.seed!(1234)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGnew;equilibration_steps=10nThermal,pre_equilibration_steps=100nThermal,scatter_fraction= 0.9) for _ in 1:20*4])
@time resultsSymm = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,psiSymm;equilibration_steps=10nThermal,pre_equilibration_steps=100nThermal,scatter_fraction= 0.9) for _ in 1:20*4])
##


# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsNaive,CT,normalize=false,dense=true,τ = 20)
# plotEnergies!(resultsOld,CT,normalize=false,dense=true,τ = 20,color = :blue)
plotEnergies!(results,CT;color=:red,nThermal = 100,normalize=false,dense=true,τ = 20) # L=15
plotEnergies!(resultsSymm,CT;color=:blue,nThermal = 100,normalize=false,dense=true,τ = 20) # L=15
plotEnergies!(resultsHighAcc,CT;color=:cyan,normalize=false,dense=true,τ = 20) # L=15

# plotEnergies!(resultsPlaq,CT;nThermal=100,p=30,color=:blue,normalize=false,dense=true,τ = 20) # L=15
# plotEnergies(results,DT.nBranch;nThermal=1,p=1000,color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
## 
current_figure()
##
plotVarEn(stochReconfRes,movavg = 10,E_exact = mean(last.(SW.getEnergies.(resultsNaive,1,50))))
##
SqsGFMCNaive = SW.getSqsGFMC(resultsNaive,1:150)
SqsGFMCSymm = SW.getSqsGFMC(resultsSymm,1:150)
SqsGFMC = SW.getSqsGFMC(results,1:150)
##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22,size = (500,400))
    ax = Axis(fig[1,1],xlabel = L"\tau",ylabel = L"\mathcal{S}^{zz}(\textbf{q})")
    
    linestyles = [:solid,:dash,:dot,:dashdot,:dashdotdot]
    colors = [:black,:red,:green,:purple,:orange]

    sqex = dropmean(SqsGFMCSymm,dims=4)[:,:,end]

    inds = sort(collect(CartesianIndices(sqex))[:], by = x -> sqex[x],rev=true)
    for (SqsGFMC,color) in zip((SqsGFMCNaive,SqsGFMCSymm),colors)

        Sqmean = dropmean(SqsGFMC,dims = 4)
        Sqerr = dropstd(SqsGFMC,dims = 4)

        # inds = [
        #     (4,5),
        #     # (5,5),
        #     # (6,5),
        #     # (5,4),
        #     # (5,6),
        #     Tuple(argmax(Sqmean[:,:,end]))
        # ]
        indsplot = inds[1:10:30]

        x = axes(SqsGFMC)[3] .* CT.τ
        for (linestyle, I) in zip(linestyles,indsplot)
            i,j = Tuple(I)
            sqm = Sqmean[i,j,:]
            sqe = Sqerr[i,j,:]
            
            l = lines!(x,sqm;linestyle,color)
            band!(x,sqm .- sqe,sqm .+ sqe ,color = (l.color[],0.2))

        end
    end
    fig
end
##
with_theme(theme_PiTicks()) do 
    Sq = dropmean(SqsGFMC,dims=4)[:,:,30] ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig = Figure(fontsize = 22,size = (800,400))
    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible = false,yticklabelsvisible=false)

    axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible = false,yticklabelsvisible=false)

    err = dropstd(SqsGFMC,dims=4)[:,:,10] ./4
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)

    fittingcoefs = optimizeCoeffs(Sq)

    SqFT = [SqFieldTheory(x,y,fittingcoefs...) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
    hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis)

    Colorbar(fig[2,1],hmMC,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,2],hmerr,label = L"\sigma( \mathcal{S}^{zz}(\textbf{q}))",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,3],hmFT,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.8),vertical=false,ticks = SimpleTicks())

    rowsize!(fig.layout,2,Relative(0.1))
    fig
end

##
let 
    Gnp
    sqtau = getImagTimeCorr(Gnp,reconfigurationTable,ObsFunc::T,mtau=size(Gnp,2)÷4, m=size(Gnp,2)÷2)
    
end