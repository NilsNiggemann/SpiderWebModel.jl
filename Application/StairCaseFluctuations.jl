import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers
import SpiderWebModel: getStairCase

##
let 
    state = getStairCase(13) 
    GC.gc()
    @time SW.getAllNeighborStates(state)
    # @time SW.generateAllPaths(state)
    
end
##
function plotPath()
    state = SW.periodicState6x6_3(12)
    @time res = SW.getAllNeighborStates(state)
    current = 1

    display(SW.plotApplPlaquettes(res.AllStates[current]))
    for i in 1:20
        current = rand(res.Neighbors[current])
        display(SW.plotApplPlaquettes(res.AllStates[current]))
    end
end
plotPath()
##

##
# Stair = getStairCase(12)
# @profview res = SW.getAllNeighborStates(Stair)
##

function solveED(state,μ)
    @time res = SW.getAllNeighborStates(state)
    @time H = SW.H(res.AllStates,res.Neighbors,μ)
    @time sol = SW.SolveH(H)
    return (;res,sol,H)
end
##
function plotMagnetization!(ax,res,sol)
    state = res.AllStates[1].parent
    mag = [SW.getMagnetization(res.AllStates,sol,CartesianIndex(i,j)) for i in axes(state,1),j in axes(state,2)]
    @info "" mavg = sum(abs,mag)/length(mag) mavgBulk = sum(abs,mag[3:end-2,3:end-2])/length(mag[3:end-2,3:end-2])

    heatmap!(ax,mag)
end

function plotStructureFac!(ax,res,sol;axkwargs=(;))
    Sq = SW.getStructureFac(SW.spinConfig.(res.AllStates),sol)

    k = Sq.Sq[1].itp.ranges[1]
    Sq_k = fetch.([Threads.@spawn Sq(kx,ky) for kx in k, ky in k])
    heatmap!(ax,k,k,Sq_k)
end

function plotOverview(res,sol;title=L"")
    with_theme(theme_SimpleTicks()) do  
        fig = Figure(;size = 0.8 .*(450,600))

        axmag = Axis(fig[1,1];title,xlabel = L"x",ylabel = L"y",aspect = 1)
        axSq = Axis(fig[2,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks(0:0.5pi:2pi),yticks = PiTicks(0:0.5pi:2pi))
        
        hm = plotMagnetization!(axmag,res,sol)
        Colorbar(fig[1,2],hm,label = L"\langle S^z \rangle")
        hm = plotStructureFac!(axSq,res,sol)
        Colorbar(fig[2,2],hm,label = L"\mathcal{S}^{zz}(\mathbf{q})")
        fig
    end

end

function getEnergy(sol)
    return sol.values[1]
end

function getEnergy(State::SW.SpinConfig)
    @time begin
        res = SW.getAllNeighborStates(State)
        μ = 0
        H = SW.H(res.AllStates,res.Neighbors,μ)
        sol = SW.SolveH(H)
        sol.values[1]
    end
    # return sol.values[1]
end

##

μ = 0
Sol5x5 = solveED(SW.periodicState5x5(14),μ)
##
plotOverview(Sol5x5.res,Sol5x5.sol,title = L"$5×5$ state, $μ = %$μ$")
##

Sol6x6 = solveED(SW.periodicState6x6(14),μ)
plotOverview(Sol6x6.res,Sol6x6.sol,title = L"$6×6$ state, $μ = %$μ$")
##
SolStair = solveED(SW.getStairCase(14),μ)
plotOverview(SolStair.res,SolStair.sol,title = L"stair state, $μ = %$μ$")
##
StairEnergy(L) = getEnergy(getStairCase(L))
# Energy5x5(L) = -L^2/10
# Energy5x5(L) = getEnergy(SW.periodicState5x5(L))
Energy5x5(L) = -length(SW.getApplicablePlaquettes(SW.periodicState5x5(L)))
Energy6x6(L) = getEnergy(SW.periodicState6x6(L))
Energy6x6_3(L) = getEnergy(SW.periodicState6x6_3(L))
Ls_small = 5:1:13
Ls_medium = 5:1:16
Ls_large = 5:1:50
StairEs = ([StairEnergy(L)/L^2 for L in Ls_small])
Ens_5x5 = ([Energy5x5(L)/L^2 for L in Ls_large])
Ens_6x6 = ([Energy6x6(L)/L^2 for L in Ls_medium])
Ens_6x6_3 = ([Energy6x6_3(L)/L^2 for L in Ls_medium])

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size = 0.8 .*(600,400))
    ax = Axis(fig[1,1],xlabel = L"L",ylabel = L"-E_0(L)/L^2",
    xscale = log10,yscale = log10
    )
    scatterlines!(ax,Ls_small,-StairEs,label = L"Staircase$$")
    scatterlines!(ax,Ls_large,-Ens_5x5,label = L"5×5")
    scatterlines!(ax,Ls_medium,-Ens_6x6,label = L"6×6")
    scatterlines!(ax,Ls_medium,-Ens_6x6_3,label = L"⋱_{6×6}")
    lines!(ax,Ls_large,ones(length(Ls_large)) ./10,label = L"E_0 =-L^2/10",linestyle = :dash, color = :grey)
    lines!(ax,Ls_large,ones(length(Ls_large)) ./16,label = L"E_0 =-L^2/16",linestyle = :dashdot, color = :grey)
    axislegend(ax,position = :rb,nbanks = 2)
    save("exactFig/EnergyScaling.png",fig)
    fig

end
