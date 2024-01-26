import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers
import SpiderWebModel: getStairCase

##


##
function plotPath()
    state = SW.periodicState5x5(12)
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
Energy5x5(L) = getEnergy(SW.periodicState5x5(L))
Energy6x6(L) = getEnergy(SW.periodicState6x6(L))
Ls = 5:1:14
StairEs = ([StairEnergy(L) for L in Ls])
Ens_5x5 = ([Energy5x5(L) for L in Ls])
Ens_6x6 = ([Energy6x6(L) for L in Ls])

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size = 0.8 .*(600,400))
    ax = Axis(fig[1,1],xlabel = L"L",ylabel = L"-E_0(L)",
    # xscale = log10,yscale = log10
    )
    scatterlines!(ax,Ls,-StairEs,label = L"Staircase$$")
    scatterlines!(ax,Ls,-Ens_5x5,label = L"5×5")
    scatterlines!(ax,Ls,-Ens_6x6,label = L"6×6")
    lines!(ax,Ls,Ls.^2 ./16,label = L"E_0 =-L²/16",linestyle = :dash, color = :grey)
    axislegend(ax,position = :rb)
    save("exactFig/EnergyScaling.png",fig)
    fig

end
