import Pkg
# Pkg.activate("../Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
meanstd(x) = (mean(x),std(x))
##

S = SW.stencilConfig(zeros(14,14),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

μ = 1.2
CT = SW.ContinuousTimeMethod(0.3,1,(1-μ)* 0.266*prod(size(S)),SW.Hxx_RK(μ))
# ψG = SW.fullVariationalFunction(S,0.15)
ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))
##
stochReconfRes = SW.stochastic_reconfiguration(S,CT,i->round(Int,2000+ 20i)÷10,ψG,200,i -> 80.,SW.IterativeSRSolver();Nwalkers = 
6*80,reconfigure=true,reset = false,rel_tolerance=0,equilibration_steps=100,pre_equilibration_steps=5_000)
plotVarEn(stochReconfRes)
##
# ψG = SW.fullVariationalFunction(S,0.1)
# ψG = SW.RKFunction()
CT = SW.ContinuousTimeMethod(0.3,1,stochReconfRes.E0[end],SW.Hxx_RK(μ))
##
ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))
@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,50,10000,ψG,equilibration_steps=2000,pre_equilibration_steps=5_000,scatter_fraction=0.5) for i in 1:12])
##

# ψG = SW.RKFunction()
ψG = SW.PlaquetteNumberGuidingFunction(only(unique(stochReconfRes.params)))
# ψG = typeof(ψG)(stochReconfRes.params)
# ψG = SW.fullVariationalFunction(S,0.1)

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,50,10000,ψG,equilibration_steps=2000,pre_equilibration_steps=5_000,scatter_fraction=0.5) for i in 1:12])
##
# plotEnergies(results,CT,p=100;normalize=true,dense=true,τ = 30,Emax = stochReconfRes.E0[end]+stochReconfRes.ΔE[end]/1.5,Emin = stochReconfRes.E0[end]-stochReconfRes.ΔE[end]/1.5)
plotEnergies(results,CT;normalize=false,dense=true,τ = 150,color = :red)
plotEnergies!(resultsOld,CT;normalize=false,dense=true,τ = 150)
current_figure()

##
mProj = round(Int,2 ÷ CT.τ)
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj) for res in results]
##
# BOp = SW.PlaquetteFlipOperator(S)
# resB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,mProj,BOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])


# BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resB,results) ])

# BValsMean = mean(BVals,dims=2)[:]
# BValsStd = std(BVals,dims=2)[:]

# scatterlines(eachindex(BValsMean).* CT.τ,BValsMean,label="plaquette")

# band!(eachindex(BValsMean).* CT.τ,BValsMean - BValsStd , BValsMean + BValsStd,color = (:black,0.2))
# current_figure()
##

RBOp = SW.RandomPlaquetteFlipOperator(S)
@time resRB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,mProj,RBOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])
##
BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resRB,results) ]) ./ length(collect(SW.plaquetteIterator(S)))

BValsMean = mean(BVals,dims=2)[:]
BValsStd = std(BVals,dims=2)[:]

scatterlines(eachindex(BValsMean).* CT.τ,BValsMean,label="plaquette",color = (:black))

band!(eachindex(BValsMean).* CT.τ,BValsMean - BValsStd , BValsMean + BValsStd,color = (:black,0.2))
current_figure()
##
Sqs = SW.getSqsGFMC(results[1:2],round(Int,10 ÷ CT.τ),CT.nBranch)
with_theme(theme_PiTicks())  do 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"q_x", ylabel = L"q_y")
    Sqm = mean(Sqs)
    Sq = SW.getSqCont(Sqm)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    hm = heatmap!(ax,qx,qy,Sq.(qs),colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = L"S(q)")
    fig
end
##
function makeBPlot(mus,energies,energiesStd,BValsMeans,BValsStds;normalizeFac = 1)
    with_theme(theme_SimpleTicks())  do 
        fig = Figure(size = (800, 600))
        axen = Axis(fig[1, 1], xlabel = L"μ", ylabel = L"E",xlabelvisible = false,xticklabelsvisible= false)
        axDE = Axis(fig[1, 1], xlabel = L"μ", ylabel = L"dE/d\mu",
            yaxisposition=:right,yticklabelcolor=:red,
            yticks = SimpleTicks(),
            xlabelvisible = false,
            ygridvisible=false,
            xgridvisible=false,
            xticklabelsvisible= false,
            xticksvisible= false,
            # xminorticksvisible=false
            ylabelvisible=true,
            ylabelcolor = :red
        )
        ax = Axis(fig[2, 1], xlabel = L"μ", ylabel = L"B")
        
        e0 = energies ./ normalizeFac
        e0Std = energiesStd ./ normalizeFac

        scatterlines!(axen,mus,e0,color = :black)
        # lines!(axen,mus,-(1 .-mus).*0.266*prod(size(S)),color = :grey,linestyle = :dash)
        linkxaxes!(axen,ax,axDE)
        perm = sortperm(mus)
        de = diff(e0[perm]) ./ diff(mus[perm])
        scatterlines!(axDE,mus[perm][1:end-1],de,color = :red)
        errorbars!(axen,mus,e0,e0Std,color = :black,whiskerwidth = 10)
        scatterlines!(ax,mus,BValsMeans,BValsStds,color = :black)
        errorbars!(ax,mus,BValsMeans,BValsStds,color = :black,whiskerwidth = 10)
        fig
    end
end
function getmupoint(S,μ,P,mProj;NSteps=1000,NStepsStochRec = NSteps,Nwalkers=20,dt=40.,NIter=50,stochRecKwargs = (;))
    ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))
    
    CT = SW.ContinuousTimeMethod(0.1,1,(1-μ)*0.266*prod(size(S)),SW.Hxx_RK(μ))

    @time stochReconfRes = SW.stochastic_reconfiguration(S,CT,NStepsStochRec,ψG,NIter,dt,SW.IterativeSRSolver();Nwalkers =60,reconfigure=true,rel_tolerance=1e-20,equilibration_steps=100,pre_equilibration_steps=1_000,verbose=false,stochRecKwargs...)
    α = only(unique(stochReconfRes.params))
    ψG = SW.PlaquetteNumberGuidingFunction(α)

    println((;μ,α))
    @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG,equilibration_steps=1000,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:6])

    # display(plotEnergies(results,CT;normalize=false,dense=true,τ = 3,color = :red))

    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj) for res in results]
    
    ens = stack(SW.getEnergies.(results,1,P))
    e_mean = reshape(mean(ens,dims=2),P)
    e_std = reshape(std(ens,dims=2),P)

    RBOp = SW.RandomPlaquetteFlipOperator(S)
    @time resRB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,mProj,RBOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])

    BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resRB,results) ]) ./ length(collect(SW.plaquetteIterator(S)))
    BValsMean = mean(BVals,dims=2)[:]
    BValsStd = std(BVals,dims=2)[:]
    
    return (;e_mean,e_std,BValsMean,BValsStd,results)
end
##
S = SW.stencilConfig(zeros(12,12),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

# mus = LinRange(-0.,0.5,4)
P = 50
mProj = round(Int,5 ÷ 0.1)
# mus = collect(1.2:-0.05:-0.1)
mus = collect(-0.1:0.05:1.2)
BValsMeans = fill(NaN,length(mus),mProj)
BValsStds = fill(NaN,length(mus),mProj)
energies = fill(NaN,length(mus),P)
energiesStd = fill(NaN,length(mus),P)


for (i,μ) in enumerate(mus)
    addkwargs = (;NSteps=1000,Nwalkers = 20)
    if μ <=0.1
        addkwargs = (;NSteps = 3000,NStepsStochRec = 3000,NIter = 100,dt = 5.,Nwalkers = 40)
    end
    (;e_mean,e_std,BValsMean,BValsStd,results) = getmupoint(S,μ,P,mProj;addkwargs...)
    energies[i,:] .= e_mean
    energiesStd[i,:] .= e_std
    BValsMeans[i,:] = BValsMean
    BValsStds[i,:] = BValsStd

    fig = makeBPlot(mus,energies[:,end],energiesStd[:,end],BValsMeans[:,end],BValsStds[:,end])
    ax = Axis(fig[3,1], ylabel = "E",xlabelvisible = false,xticklabelsvisible= false)
    ax2 = Axis(fig[4,1],xlabel = L"projection $$")
    linkxaxes!(ax,ax2)
    CT = SW.ContinuousTimeMethod(0.1,1,(1-μ)*0.266*prod(size(S)),SW.Hxx_RK(μ))

    plotEnergies!(ax,results,CT,p=100;normalize=true,dense=true,τ = 8,color = :black)
    lines!(ax2,(axes(BValsMean,1) .-1) .*CT.τ,BValsMean,color = :red)
    band!(ax2,(axes(BValsMean,1) .-1) .*CT.τ,BValsMean .- BValsStd,BValsMean .+ BValsStd ,color = (:red,0.2))
    display(fig)
end
##
makeBPlot(mus[2:end],energies[2:end,10],energiesStd[2:end,10],BValsMeans[2:end,10],BValsStds[2:end,10];normalizeFac = length(S))
##
# for (i,μ) in enumerate(mus[1:6])
for μ in [0.25]
    i = findfirst(isequal(μ),mus)
    (;e_mean,e_std,BValsMean,BValsStd,results) = getmupoint(S,μ,P,mProj,NSteps=1000,Nwalkers = 120,dt = 20,NStepsStochRec = i->(200 +1i),NIter = 100,stochRecKwargs = (;Nwalkers = 120))
    # (;e_mean,e_std,BValsMean,BValsStd,results) = getmupoint(S,μ,P,mProj,NSteps=1000,Nwalkers = 20,dt = 0.4,NStepsStochRec = 3000,NIter = 100)
    energies[i,:] .= e_mean
    energiesStd[i,:] .= e_std
    BValsMeans[i,:] = BValsMean
    BValsStds[i,:] = BValsStd


    fig = makeBPlot(mus,energies[:,end],energiesStd[:,end],BValsMeans[:,end],BValsStds[:,end])
    ax = Axis(fig[3,1], ylabel = "E",xlabelvisible = false,xticklabelsvisible= false)
    ax2 = Axis(fig[4,1],xlabel = L"projection $$")
    linkxaxes!(ax,ax2)
    CT = SW.ContinuousTimeMethod(0.1,1,(1-μ)*0.266*prod(size(S)),SW.Hxx_RK(μ))

    plotEnergies!(ax,results,CT,p=100;normalize=true,dense=true,τ = 8,color = :black)
    lines!(ax2,(axes(BValsMean,1) .-1) .*CT.τ,BValsMean,color = :red)
    band!(ax2,(axes(BValsMean,1) .-1) .*CT.τ,BValsMean .- BValsStd,BValsMean .+ BValsStd ,color = (:red,0.2))
    display(fig)
end