import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
##
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# S = SW.stencilConfig(SW.constructConfigPath(15,15,SW.ALLGS_S12),1/2)
HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1][:,1]
HConfs = SW.spinConfig.(HStair.AllStates,Ref(S),Ref(HStair.plaqMapping))
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf)
        ind = AllSTDict[Conf]
        return v0[ind] + 0.003
    end
end
##


nThermal = 400_000
results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1500_000,SW.ConstructVaritationalFunc(0.197,S),3) for _ in 1:8])
# results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+500_000,constructExactGuidingFunc(v0,HConfs),3) for _ in 1:8])
##
ens = [SW.getEnergies(res.TotalWeights,res.energies,nThermal,150) for res in results]
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
    enmat = stack(ens)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax)
    xlims!(ax,1,length(en))
    ylims!(ax,E0-1e-2,E0+1e-1)
    fig
end
##