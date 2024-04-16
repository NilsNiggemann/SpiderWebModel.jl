import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
##
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# S = SW.stencilConfig(SW.constructConfigPath(15,15,SW.ALLGS_S12),1/2)
HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1][:,1]
HConfs = SW.spinConfig.(HStair.AllStates,Ref(SW.
SpinConfig(S)),Ref(HStair.plaqMapping))
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf)
        ind = AllSTDict[Conf]
        return v0[ind] + 0.01
    end
end
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
function binningDataEstimates(resultsVec,binsize)

    energies = Vector{Vector{Float64}}()
    TotalWeights = Vector{Vector{Float64}}()

    for res in resultsVec
        Enbins = SW.splitIntoBins(res.energies,binsize)
        bxbins = SW.splitIntoBins(res.TotalWeights,binsize)
        append!(energies,Enbins)
        append!(TotalWeights,bxbins)
    end
    return (;energies,TotalWeights)
end
# function binningDataEstimates(resultsVec,binsize)

#     newResVec = empty(resultsVec[1])

#     for res in resultsVec
#         bins = SW.splitIntoBins(res,binsize)
#         bxbins = SW.splitIntoBins(res.TotalWeights,binsize)
#         append!(energies,Enbins)
#         append!(TotalWeights,bxbins)
#     end
#     return (;energies,TotalWeights)
# end
##
nThermal = 4_000
# SW.Random.seed!(12345)
# results = [SW.startManyWalkerGFMC(S,2,55_000,3,nThermal,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:35]
results = [SW.startManyWalkerGFMC(S,4,1_000_000÷4,4,nThermal,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:12]
# @time results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1000_000,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:24])

##
# @time obs = [SW.getObservables(res,S,float,1,100) for res in results]
ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,100) for res in results]
# binnedData = binningDataEstimates(results,100nThermal)
# @time ens = fetch.([Threads.@spawn SW.getEnergies(binnedData.TotalWeights[i],binnedData.energies[i],1,100) for i in eachindex(binnedData.TotalWeights)])
# locEn = collect(SW.splitIntoBins(results[1].energies,100000))
# locW = collect(SW.splitIntoBins(results[1].TotalWeights,100000))
# @time ens = fetch.([Threads.@spawn SW.getEnergies(locW[i],locEn[i],1,100) for i in eachindex(locEn)])
en = mean(ens)
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    # ens = getfield.(obs,:E0)
    # en = mean(ens)
    err = sqrt.(var(ens))
    scatter!(ax,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,eachindex(en),en,err,whiskerwidth = 3.5,color = :black)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax,merge=true)
    xlims!(ax,0.5,length(en))
    ylims!(ax,E0-1e-2,E0+1e-1)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end