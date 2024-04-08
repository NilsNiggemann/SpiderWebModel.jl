import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
##
HStair = SW.generateHilbertSpace(SW.getStairCase(12))
E0 = SW.SolveHKrylov(HStair.H).values[1]
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
##
nThermal = 400_000
results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1_200_000,SW.ConstructVaritationalFunc(0.197,S),2) for _ in 1:16])
##
en = mean([SW.getEnergies(res.TotalWeights,res.energies,nThermal,150) for res in results])
# en = mean(
#     fetch.([Threads.@spawn [SW.getEnergy(res.TotalWeights,res.energies,p,20_000) for p in 1:80] for res in results])
# )
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0")
    scatter!(ax,en,label = L"GFMC$$",color = :black, marker = '×',markersize = 14)
    hlines!([E0],color = :red,label = L"exact $$")
    axislegend(ax)
    fig
end
##
function test()
    S = SW.stencilConfig(parent(SW.getStairCase(9)),1/2)
    mv = filter(!=((0,0,0)),results[1].Allmoves)
    allSt = Set([collect(S)])
    for m in mv
        SW.applyPlaquette!(S,m...)
        if S ∉ allSt
            push!(allSt,collect(S))
        end
    end
    allSt
end
test()
##
# ##
# for i in HStair.AllStates
#     S = SW.spinConfig(i,SW.getStairCase(8),HStair.plaqMapping)
#     SW.plotApplPlaquettes(S) |> display
# end