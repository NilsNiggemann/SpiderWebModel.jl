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

S = SW.stencilConfig(zeros(18,18),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
# S .= SW.h5read("temp.h5","conf")
μ = 0.8
# ψG = SW.fullVariationalFunction(S,0.15*(1-μ))
ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))

##
CTFindOpt = SW.ContinuousTimeMethod(10.,1,(1-μ)* 0.2*0.266*length(S),SW.Hxx_RK(μ))

@time OptimStart = SW.startManyWalkerGFMC(S,CTFindOpt,30,100,ψG;equilibration_steps=1,pre_equilibration_steps=50_000,scatter_fraction=0.5)
OptimStart.TotalWeights
ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))
##
initializer = SW.WeightedConfigsInitializers(OptimStart.SaveConfigs,OptimStart.TotalWeights)
CT_SR = SW.ContinuousTimeMethod(0.05,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))

stochReconfRes = SW.stochastic_reconfiguration(S,CT_SR,i->round(Int,100+ 2i),ψG,50,i->100*min(0.4,0.05+0.001*i),SW.IterativeSRSolver();Nwalkers = 
60,reconfigure=true,reset = true,rel_tolerance=0,equilibration_steps=0,initializer,
report_steps = 2,
)
plotVarEn(stochReconfRes)
##
# ψG = SW.fullVariationalFunction(S,0.1)
# ψG = SW.RKFunction()
CT = SW.ContinuousTimeMethod(0.1,1,stochReconfRes.E0[end],SW.Hxx_RK(μ))
##
ψG = SW.localPlaquetteGuidingFunction(S,0.15*(1-μ))
# ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))

@time resultsOld = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,80,3000,ψG,equilibration_steps=1000,pre_equilibration_steps=10_000,scatter_fraction=0.0) for i in 1:12])
##

# ψG = SW.RKFunction()
# ψG = SW.PlaquetteNumberGuidingFunction(only(unique(stochReconfRes.params)))
# ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
# ψG = typeof(ψG)(stochReconfRes.params)
ψG = SW.LocalPlaquetteGuidingFunction(stochReconfRes.params)
# ψG = SW.FullVariationalGuidingFunction(stochReconfRes.params)
# ψG = SW.fullVariationalFunction(S,0.1)
CT2 = SW.ContinuousTimeMethod(0.1,1,-stochReconfRes.E0[end],SW.Hxx_RK(μ))

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,80,3000,ψG;equilibration_steps=1000,initializer) for i in 1:12])
##
# plotEnergies(results,CT,p=100;normalize=true,dense=true,τ = 30,Emax = stochReconfRes.E0[end]+stochReconfRes.ΔE[end]/1.5,Emin = stochReconfRes.E0[end]-stochReconfRes.ΔE[end]/1.5)
plotEnergies(results,CT2;normalize=false,dense=true,τ = 30,color = :red)
plotEnergies!(resultsOld,CT;normalize=false,dense=true,τ = 30)
current_figure()
##
equilib_plots(results;scatter_fraction=0,averageSteps=10,Ntrack=30,p = round(Int,20÷CT.τ),plotPopulation=true)

##
SqsGFMC = fetch.([Threads.@spawn SW.getSqGFMC(res,round(Int,50÷CT.τ)+2) for res in results])

function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end

function SqFieldTheory2(kx,ky,k,b2)
    (sqrt(2)*sqrt(-((-4 + 4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky))*(40*b2 + k + 8*b2*(cos(2*kx) + 2*cos(kx)*(-4*cos(ky) + cos(kx)*cos(2*ky))))))*   (cos(kx) - cos(ky) + 2*sin(kx)*sin(ky))^2)/((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2+1e-30)
end

function makeSqFTPlots(SqsGFMC,k=1,b2=0)
    SqFT_func(x,y)  = SqFieldTheory2(x,y,k,b2)
    with_theme(theme_PiTicks()) do 
        # Sq = sqrt.(var(real(SqsGFMC))) ./4
        Sq = mean(SqsGFMC) ./4
        kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
        fig = Figure(fontsize = 22,size = (800,400))
        axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
        axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        # axDiff = Axis(fig[1,4],xlabel = L"k_x",ylabel = L"k_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

        # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
        err = sqrt.(var(SqsGFMC)) ./4
        hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
        SqFT = [SqFT_func(x,y) for x in kx, y in ky]
        hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
        # heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(!isnan,Sq))
        hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis)
        # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)

        Colorbar(fig[2,1],hmMC,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
        Colorbar(fig[2,2],hmerr,label = L"\sigma( \mathcal{S}^{zz}(\textbf{q}))",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
        Colorbar(fig[2,3],hmFT,label = L" \mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.8),vertical=false,ticks = SimpleTicks())

        rowsize!(fig.layout,2,Relative(0.1))
        fig
    end
end
makeSqFTPlots(SqsGFMC,1,0.1)

##
function plotReconfStats(reconfigurationTable)
    survWalkers = SW.StatsBase.countmap(length.(unique.(eachcol(reconfigurationTable))))
    surv = sort(collect(keys(survWalkers)))
    counts = [survWalkers[s] for s in surv] ./ size(reconfigurationTable,2)
    barplot(surv,counts;axis = (;xlabel = "Surviving walkers",ylabel = "Count"),color = :black)
end
plotReconfStats(results[3].reconfigurationTable[:,1000:end])
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
Sqs = SW.getSqsGFMC(results,round(Int,5 ÷ CT.τ),CT.nBranch)
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