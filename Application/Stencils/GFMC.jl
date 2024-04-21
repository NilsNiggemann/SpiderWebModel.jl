import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
##
S = SW.stencilConfig(parent(SW.getStairCase(8)),1/2)
# S = SW.stencilConfig(SW.constructConfigPath(15,15,SW.ALLGS_S12),1/2)
HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1]
HConfs = SW.spinConfig.(HStair.AllStates,Ref(SW.
SpinConfig(S)),Ref(HStair.plaqMapping))
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf)
        ind = AllSTDict[Conf]
        return v0[ind] + 0.01
    end
end

function makeBlocks(arr;blocksize = nothing, numBlocks = 32)
    if blocksize === nothing
        blocksize = length(arr) ÷ numBlocks
    end
    blockedArr = collect(SW.splitIntoBins(arr,blocksize))
    
    minsize,maxsize = extrema(length.(blockedArr))
    
    if maxsize!=minsize
        pop!(blockedArr)
    end
    return blockedArr
end
magEx = SW.getMagnetization(HConfs, v0)
##


nThermal = 4_000
results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+250_000,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:12*12])
# results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1300_000,SW.ConstructVaritationalFunc(0.197,S),3) for _ in 1:4*4])
# results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+900_000,constructExactGuidingFunc(v0,HConfs),3) for _ in 1:8])
##
ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,100) for res in results]
en = mean(ens)

# en = SW.getEnergies(results[6].TotalWeights,results[6].energies,nThermal,150)
# en = mean(
#     fetch.([Threads.@spawn [SW.getEnergy(res.TotalWeights,res.energies,p,20_000) for p in 1:80] for res in results])
# )
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    scatter!(ax,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,eachindex(en),en,sqrt.(var(ens)),whiskerwidth = 3.5,color = :black)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax)
    xlims!(ax,0.5,length(en))
    ylims!(ax,E0-1e-2,E0+1e-2)
    fig
end
##

@time obs = fetch.([Threads.@spawn SW.getObservables(res,S,float,nThermal,100) for res in results])
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    en = mean(getfield.(obs,:E0))
    err = sqrt.(var(getfield.(obs,:E0)))
    scatter!(ax,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,eachindex(en),en,err,whiskerwidth = 3.5,color = :black)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax,merge=true)
    xlims!(ax,1,length(en))
    ylims!(ax,E0-1e-2,E0+1e-1)
    fig
end
#
with_theme(theme_SimpleTicks()) do
    # fig,ax,hm = heatmap(sqrt.(var(getfield.(obs,:Obs) ./ 2)),colormap = :grays,axis=(;aspect=1))
    fig,ax,hm = heatmap(mean(getfield.(obs,:Obs) ./ 2),colormap = :grays,axis=(;aspect=1))
    # Colorbar(fig[1,2],hm)
    # Colorbar(fig[1,2],hm,ticks = ([-0.5,0.,0.5],[L"|\downarrow>",L"0",L"|\uparrow>"]))
    Colorbar(fig[1,2],hm,ticks = ([-0.5,0.,0.5]))
    fig
end
#
#___________ManyWalkers_______________________
##

nThermal = 2_000
SW.Random.seed!(1234)
# results = [SW.startManyWalkerGFMC(S,2,55_000,3,nThermal,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:35]
nBra = 3
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,16,nThermal+400_000÷nBra,nBra,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:32])
##
# nBra = 1
# @time results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1200_000,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:6*4])
# ens = [SW.getEnergies(res.TotalWeights[nThermal:end],res.energies[nThermal:end],1,150÷nBra) for res in results]
ens = [SW.getEnergies(w,e,1,150÷nBra) for res in results[1:end] for (w,e) in zip(makeBlocks(res.TotalWeights[nThermal:end],numBlocks=2),makeBlocks(res.energies[nThermal:end],numBlocks=2)) ]

en = mean(ens)
# errs = getErrBlocking(results[1].energies[nThermal:end],results[1].TotalWeights[nThermal:end],2*10^4,20,E0) ./ ( (length(results[1].energies)-nThermal) ÷ 2*10^4)
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    # ens = getfield.(obs,:E0)
    en = mean(ens)
    # M = length(results[1].energies)
    # Mk = M ÷ length(ens)
    # println(Mk)
    err = sqrt.(var(ens))
    # err = 0.004 .* ones(length(en))
    proj = nBra .*eachindex(en)
    scatter!(ax,proj,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,proj,en,err,whiskerwidth = 3.5,color = :black)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax,merge=true)
    xlims!(ax,0.5,last(proj))
    ylims!(ax,E0-1e-2,E0+3e-2)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end

##
  # @time obs = [SW.getObservables(res,S,float,1,150 ÷nBra) for res in results]
# ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,150÷nBra) for res in results]
# binnedData = binningDataEstimates(results,100nThermal)
# @time ens = fetch.([Threads.@spawn SW.getEnergies(binnedData.TotalWeights[i],binnedData.energies[i],1,100) for i in eachindex(binnedData.TotalWeights)])
function getEnsBinning(energies,weights,nThermal,BinSize)
    @views locEn = collect(SW.splitIntoBins(energies[nThermal:end],BinSize))
    @views locW = collect(SW.splitIntoBins(weights[nThermal:end],BinSize))
    # for res in results[2:end]
    #     append!(locEn,collect(SW.splitIntoBins(res.energies[nThermal:end],BinSize)))
    #     append!(locW,collect(SW.splitIntoBins(res.TotalWeights[nThermal:end],BinSize)))
    # end
    pop!(locEn)
    pop!(locW)
    @time ens = fetch.([Threads.@spawn SW.getEnergies(locW[i],locEn[i],1,150÷nBra) for i in eachindex(locEn)])
end
##
ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,150÷nBra) for res in results]
binsize = 500000

##
function getErrBlocking(energies,weights,blocksize,LProjection,E0)

    locEn = collect(SW.splitIntoBins(energies,blocksize))
    locW = collect(SW.splitIntoBins(weights,blocksize))
    
    minsize,maxsize = extrema(length.(locEn))
    
    if maxsize!=minsize
        pop!(locEn)
        pop!(locW)
    end
    # @info "" unique(length.(locEn))
    M_b = length(locEn)

    X²s = [last(SW.getEnergies(locW[i],locEn[i].^2,1,LProjection)) for i in 1:M_b]
    # safesum(x) = x<0 ? NaN : x
    # Xj = [last(SW.getEnergies(locW[i],locEn[i],1,LProjection)) for i in 1:M_b]
    σb = [√(1/(M_b-1) * sum(X² .- E0^2)) for X² in X²s]
    # σb = [√(1/(M_b-1) * sum((X .- E0)^2)) for X in Xj]


    # err = sqrt(var(last_en))
 
    # err = sqrt()
end

with_theme(theme_SimpleTicks()) do
    # BlockSizes = [2^i for i in 3:10]
    maxSize = length(results[1].energies) -nThermal
    BlockSizes = [round(Int,2^i) for i in LinRange(log2(100),log2(maxSize ÷4),40)]
    E0MC = ens[1][20÷nBra]
    
    err = fetch.([Threads.@spawn last(getErrBlocking(results[1].energies[nThermal:end],results[1].TotalWeights[nThermal:end],bs,20÷nBra,E0MC)) for bs in BlockSizes])
    
    numBins = maxSize .÷ BlockSizes
    # return numBins
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"block size $$",ylabel = L"\sigma_b",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5),xscale = log10)

    # RX = BlockSizes .* err.^2 ./1
    # scatter!(ax,BlockSizes,RX)
    scatter!(ax,BlockSizes,err)
    vlines!(ax,BlockSizes[findfirst(<(32),numBins)],color = :red)
    fig
    # scatter(BlockSizes,E_mean)
    # hlines!(E0)
    # current_figure()
end
# blockingAnalysis(results[1].energies[nThermal:end],results[1].TotalWeights[nThermal:end],2^7-1,150÷nBra)
##
using BinningAnalysis

function getEJackKnife(TotWeight,energies,LProjection)
    Gnp = SW.precomputeNormalizedAccWeight(TotWeight,1,LProjection)
    
    # inner(energies,weights) = SW.getEnergies(weights,energies,1,LProjection;Gnp)

    # inner(weights,energies) = SW.getEnergies(weights,energies,1,LProjection)
    # e,emean = jackknife(inner,TotWeight,energies)
end
getEJackKnife(results[1].TotalWeights[nThermal:end],results[1].energies[nThermal:end],20÷nBra)
##
# ens = getEnsBinning(results[1].energies,results[1].TotalWeights,nThermal,2*10^4)
#___________Observables_______________________
##
# getAvgMag(Conf) = sum(Conf) ./ (2*length(Conf))
function getMag(results,pmax,I)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,pmax) for res in results]

    getMag(Conf) = Conf[I] /2

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs[nThermal:end],res.reconfTable[:,nThermal:end],getMag,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
##
mags1 = getMag(results,100÷nBra,CartesianIndex(3,3)) 
mags2 = getMag(results,100÷nBra,CartesianIndex(4,1)) 
mags3 = getMag(results,100÷nBra,CartesianIndex(2,3)) 

##

function errorBarLegend(size = 0.5;linekwargs = (;),markerkwargs = (;))
    center = 0.5
    ymin = center - size/2
    ymax = center + size/2

    errorbar = [Point2f(center, ymin), Point2f(center, ymax)]
    [
        LineElement(linepoints = errorbar),
        MarkerElement(points = errorbar, marker = :hline, markersize = 10;linekwargs...),
        MarkerElement(points = [Point2f(center, center),], marker = '●', markersize = 7;markerkwargs...)
    ]
end

with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"magnetization $$",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    proj = nBra .*eachindex(mags)
    
    Cols = (:blue,:green,:orange)
    Is = ((3,3),(4,1),(2,3))
    for (i,mags) in enumerate((mags1,mags2,mags3))
        mExact = magEx[Is[i]...]

        sgn = sign(mExact)
        m =  sgn .* mean.(mags)
        scatter!(ax,proj,m,label = L"GFMC$$",color = Cols[i], marker = '●',markersize = 5)
        err =sqrt.(var.(mags))
        errorbars!(ax,proj,m,err,whiskerwidth = 3.5,color = Cols[i])
        
        hlines!([sgn * mExact],color = Cols[i],label = L"exact $$")
        sgstring = sgn > 0 ? "" : "-"
        text!(ax,(proj[end÷2],mExact-0.015),color = Cols[i],text = L"%$sgstring \langle S_{%$(Is[i])}\rangle")
    end

    Legend(fig[1, 1], [[
        errorBarLegend(0.6)
    ],
    [LineElement(linepoints = [Point2f(0., 0.5), Point2f(1, 0.5)])]], [L"GFMC$$",L"exact $$"], tellheight = false, tellwidth = false,halign = :right, valign = :center,margin = (10,10,200,10))


    # axislegend(ax,merge=true,position = :rc,unique=true)
    # xlims!(ax,0.5,last(proj))
    # ylims!(ax,E0-1e-2,E0+3e-2)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end