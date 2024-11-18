import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
# using MKL
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(12,12),1;
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
ψG = SW.JastrowFunction(S,Float64)
# ψG = SW.orderGuidingFunction(S)
# ψG = SW.PlaquetteRBM(S,1,Float64)
Symmetry = SW.SymmetryGroup(SW.TranslationalSymmetry(SA[2,2],SA[-2,2]),SW.ExchangeSymmetry())
# Symmetry = SW.SymmetryGroup(SW.ExchangeSymmetry())
ψGSymm = SW.getNonSymmetric(ψG)
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
nThermal = 1000
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
# DT = SW.DiscreteTimeMethod(0.,3,0.266*length(S))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.2))
CTSR = SW.ContinuousTimeMethod(8,Hxx = CT.Hxx)


##
cappedGrowth(x,start,stop,offset,growth) = start + 0.5(stop-start)* (1 +tanh(growth *(x -offset)))
SRSteps = 500

numSteps(i) = round(Int,cappedGrowth(i,20,40,SRSteps - SRSteps÷2,10/SRSteps))
learningRate(i) = cappedGrowth(i,8e-1,4e-1,SRSteps - SRSteps÷2,5/SRSteps)
lines(1:SRSteps,5000*learningRate.(1:SRSteps))
lines!(1:SRSteps,numSteps.(1:SRSteps))
# current_figure()
##
# learningRate(i) = i > 450 ? 1e-3 : 1e-3
# numSteps(i) = i > SRSteps ÷2 ? 30 : 20
##

stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,numSteps ,ψGSymm,SRSteps,8e-4,SW.IterativeSRSolver();Nwalkers = 1*28,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 10,
reset = false,
# outfile = "tempSR/SR2.h5"
)

ψGnew = deepcopy(ψG)
SW.get_params(ψGnew) .= stochReconfRes.params
# SW.get_params(ψGnew) .= stochReconfRes.params_steps[:,begin]
plotVarEn(stochReconfRes,movavg = 10)
##
with_theme(theme_SimpleTicks()) do
    S = SW.stencilConfig(zeros(12,12),1;
    boundaryCondition = :periodic
    )
    ψG = SW.JastrowFunction(S)

    params = SW.get_params(ψG)

    vij = SW.PeriodicMatrix(reshape(ψGnew.v_ij[1+2*size(S,1)+2,:],size(S)))
    # vij = SW.PeriodicMatrix(reshape(ψGnew.v_ij[1+2*size(S,1)+2,:],size(S)))[3:end+2,3:end+2]

    fig,ax,hm = SW.heatmap(vij)
    # use text! to annotate the heatmap with the value of each vij matrix element, recast as an Int
    # for J in CartesianIndices(vij)
    #     for I in CartesianIndices(vij)
    #         if vij[I] == vij[J] && I != J
    #             scatter!(ax, [Point(Tuple(I))])
    #         end
    #     end
    # end
    
    for I in CartesianIndices(vij)
        v = vij[I]
        i,j = Tuple(I)

        txtcolor = vij[i, j] < mean(vij) ? :white : :black

        # val = Int(vij[i,j])
        val = (vij[i,j])
        # text!(ax, "$(val)", position = (i, j),
            # color = txtcolor, align = (:center, :center))

        # text!(ax, [Point(i,j)], string(round(Int, v)), color = :white,fontsize = 25)
    end
    fig
end

##

SW.Random.seed!(1234)
NWalkers = 28
NSteps = 2000
@time resultsNaive = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGold;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])

# @time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψG;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])
##
SW.Random.seed!(1234)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])
##


# plotEnergies(results,nBra,-20.35;Emin=-20.5,Emax=-19.8) # L=10
plotEnergies(resultsNaive,CT,normalize=false,dense=true,τ = 20)
# plotEnergies!(resultsOld,CT,normalize=false,dense=true,τ = 20,color = :blue)
plotEnergies!(results,CT;color=:red,nThermal = 100,normalize=false,dense=true,τ = 20) # L=15
# plotEnergies!(resultsPlaq,CT;nThermal=100,p=30,color=:blue,normalize=false,dense=true,τ = 20) # L=15
# plotEnergies(results,DT.nBranch;nThermal=1,p=1000,color=:red) # L=15
current_figure()
# plotEnergies(results,nBra,-49.7;Emin=-50.5,Emax=-46)
## 
@time resultsHighAcc = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,5NWalkers,2NSteps,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:28])
plotEnergies!(resultsHighAcc,CT;color=:darkblue,normalize=true,dense=true,τ = 20) # L=15
current_figure()
##
plotVarEn(stochReconfRes,movavg = 10,E_exact = mean(last.(SW.getEnergies.(resultsNaive,1,50))))
##
SqsGFMCNaive = SW.getSqsGFMC(resultsNaive,1:50)
SqsGFMC = SW.getSqsGFMC(results,1:50)
##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22,size = (500,400))
    ax = Axis(fig[1,1],xlabel = L"\tau",ylabel = L"\mathcal{S}^{zz}(\textbf{q})")
    
    linestyles = [:solid,:dash,:dot,:dashdot,:dashdotdot]
    colors = [:black,:red,:green,:purple,:orange]
    for (linestyle,SqsGFMC,color) in zip(linestyles,(SqsGFMC,SqsGFMCNaive),colors)

        Sqmean = dropmean(SqsGFMC,dims = 4)
        Sqerr = dropstd(SqsGFMC,dims = 4)

        inds = [
            (4,5),
            # (5,5),
            # (6,5),
            # (5,4),
            # (5,6),
            Tuple(argmax(Sqmean[:,:,end]))
        ]

        x = axes(SqsGFMC)[3] .* CT.τ
        for (i,j) in inds
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
    Sq = dropmean(SqsGFMC,dims=4)[:,:,10] ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig = Figure(fontsize = 22,size = (800,400))
    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

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