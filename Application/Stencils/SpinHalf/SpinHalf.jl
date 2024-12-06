import Pkg
cd(@__DIR__)
Pkg.activate("../../")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("../plottingUtils.jl")
##
#___________Periodic Boundaries_______________________
function getSectorConfig(L,i)
    S = SW.stencilConfig(0.5*ones(L,L),1/2;boundaryCondition = :periodic)

    S.= SW.getSelectedS12PeriodicState(L,i)
end

##
function initializeGWF(S,Symmetry) 
    ψG = SW.SimpleJastrowFunction(S)
    ψGSymm = SW.symmetrize(ψG,Symmetry,S)
    SW.rand!(ψGSymm,1e-5)
    return ψGSymm
end

function optimizeWF!(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,NConfs,ψG,NSteps,dt,SW.IterativeSRSolver();Nwalkers = 20,reconfigure = false,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=25,kwargs...)
    SW.get_params(ψG) .= stochReconfRes.params
    if showplot
        display(plotVarEn(stochReconfRes))
    end
    return ψG
end

##

function getEnergies(S,mus,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,dt = 5e-3,NConfs = 20,nopt1 = 100,nopt2 = 30)
    ψG = initializeGWF(S,Symmetry)
    ens = zeros(length(mus),NRuns)
    for (i,mu) in enumerate(mus)
        println("mu = $mu")
        CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(mu))
        CTSR = SW.ContinuousTimeMethod(10*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

        nopt = i == 1 ? nopt1 : nopt2
        @time optimizeWF!(ψG,S,CTSR,nThermal;NSteps = nopt,NConfs,dt)
        @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG.psi) for _ in 1:NRuns])
        # display(plotEnergies(results,CT))
        en = stack(SW.getEnergies.(results,1,50))[end,:]
        ens[i,:] .= en

    end
    return ens
end

function mainRun(Ls,muRange,Symmetries)
    run(i) = getEnergies( getSectorConfig(Ls[i],i),muRange,Symmetries[i])

    allE = fetch.([Threads.@spawn run(i) for i in 1:3])

    return allE
end
##

muRange = reverse(LinRange(0,0.99,20))

Ls = [20,24,20,18]
Symms = [
    SW.TranslationalSymmetry([-2,2],[2,2]),
    SW.TranslationalSymmetry([6,0],[0,6]),
    SW.TranslationalSymmetry([5,1],[0,4]),
    SW.TranslationalSymmetry([3,-1],[0,6])
]
ens1 = getEnergies( getSectorConfig(Ls[1],1),muRange,Symms[1],NConfs=20,dt = 2e-4)
##
ens2 = getEnergies( getSectorConfig(Ls[2],2),muRange,Symms[2],NConfs=20,dt = 2e-4)
##
ens3 = getEnergies( getSectorConfig(Ls[3],3),muRange,Symms[3],NConfs=20,dt = 2e-4)
##
ens4 = getEnergies( getSectorConfig(Ls[4],4),muRange,Symms[4],NConfs=20,dt = 2e-4)
# ens2 = mainRun([20,20,20,20],muRange)
##
ens = [ens1,ens2,ens3,ens4]
E_trivial(mu) = -0.1*(1-mu)

with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"\mu",ylabel = L"E/L^2")
    axkwargs = SW.getConfigAxis(getSectorConfig(8,1))


    colors = [:blue,:green,:red,:purple]
    spincolors(color) = (topspinecolor = color,bottomspinecolor = color,leftspinecolor = color,rightspinecolor = color)

    inax = [
        insetAtPoint(fig,ax,(0.05 +0.18(i-1),-0.01),(36,36);
        spincolors(colors[i])...,
        spinewidth = 4,
        xticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        axkwargs...
        ) for i in eachindex(ens)
    ]

    inaxTriv = insetAtPoint(fig,ax,(0.9,-0.09),(50,50);
        spincolors(:grey)...,
        spinewidth = 4,
        xticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        title = L"Trivial $$",
        axkwargs...
    )

    SW.plotApplPlaquettes!(inaxTriv,getSectorConfig(20,5),markersize = 8)

    etriv = E_trivial.(muRange)# ./ (1 .-muRange)

    lines!(ax,muRange,etriv,color = :grey,linewidth = 3)
    # ylims!(ax,-0.12,0.03)
    for (i,en) in enumerate(ens)
        L = Ls[i]
        SW.plotApplPlaquettes!(inax[i],getSectorConfig(L,i),markersize = 8)
        color = colors[i]
        e = dropmean(en,dims=2) ./ L^2 # ./ (1 .-muRange)
        e_err = dropstd(en,dims=2) ./ L^2 # ./ (1 .-muRange)
        scatterlines!(ax,muRange,e;color,marker = '×',markersize = 15)
        errorbars!(ax,muRange,e,e_err;color)
        # errlines!(ax,muRange,e,e_err;color)
    end
    fig
end

##

S = getSectorConfig(24,1)
ψG = initializeGWF(S,SW.TranslationalSymmetry([-2,2],[2,2]))
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(0.0))
CTSR = SW.ContinuousTimeMethod(30*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

@time optimizeWF!(ψG,S,CTSR,1000;NSteps = 500,NConfs=20,dt=2e-4,showplot=true,report_steps = 2)

##
results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20*4,2000,ψG.psi) for _ in 1:10])

Sqs = SW.getSqsGFMC(results,1:100)

plotEnergies(results,CT)
##


with_theme(theme_PiTicks()) do
    fig = Figure(fontsize = 22,size = (400,320))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks((0,pi)),yticks = PiTicks((0,pi)))
    # ax2 = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    SqMat = dropmean(Sqs,dims=4)[:,:,end]
    SqError = dropstd(Sqs,dims=4)[:,:,end]

    # fittingCoefs = optimizeCoeffs(SqMat)
    Sq = SW.getSqCont(SqMat)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    hm = heatmap!(ax,qx,qy, Sq.(qs))
    # hm = heatmap!(ax,qx,qy, Sq.(qs) .+1e-10,colorscale = log10,colorrange = (0.01,maximum(SqMat)))
    Colorbar(fig[1,2],hm,label = L"\langle S^z(\mathbf{q})^*  S^z(\mathbf{q})\rangle")
    # SQFT(x) = SqFieldTheory(SVector(x),fittingCoefs)
    # heatmap!(ax2,qx,qy,SQFT.(qs))
    fig
end

##
S = getSectorConfig(24,1)
ψG = SW.RKFunction()

CT = SW.ContinuousTimeMethod(0.5,w_avg_estimate = 0.,Hxx = SW.Hxx_RK(0.0))

resultsRK = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20*1,8000,ψG,pre_equilibration_steps = 100000) for _ in 1:20]) 
##
SqsRK = SW.getSqsGFMC(resultsRK,1:100)
with_theme(theme_PiTicks()) do
    fig = Figure(fontsize = 22,size = (400,320))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks((0,pi)),yticks = PiTicks((0,pi)))
    SqMat = dropmean(SqsRK,dims=4)[:,:,end]
    SqError = dropstd(SqsRK,dims=4)[:,:,end]

    Sq = SW.getSqCont(SqMat)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    hm = heatmap!(ax,qx,qy, Sq.(qs))
    Colorbar(fig[1,2],hm,label = L"\langle S^z(\mathbf{q})^*  S^z(\mathbf{q})\rangle")
    fig
end