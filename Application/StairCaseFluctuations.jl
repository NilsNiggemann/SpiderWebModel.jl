import SpiderWebModel as SW
using StaticArrays, CairoMakie
using MakieHelpers
import SpiderWebModel: getStairCase
##
function getPeriodic(parent)
    state = parent |> Array
    SW.SpinConfig(SW.PeriodicMatrix(state),parent.S)
end
##
let 
    state = getStairCase(13) 
    GC.gc()
    @time SW.getAllNeighborStates(state)
    # @time SW.generateAllPaths(state)
    
end
##
function plotPath()
    # state = SW.getStairCase(8)
    state = getPeriodic(SW.periodicState5x5(10))
    # state = SW.SpinConfig(SW.PeriodicMatrix(state),1/2)

    @time res = SW.getAllNeighborStates(state)
    current = 1

    display(SW.plotApplPlaquettes(res.AllStates[current]))
    for i in 1:20
        current = rand(res.Neighbors[current])
        display(SW.plotApplPlaquettes(res.AllStates[current]))
    end
    state
end
plotPath()
##

##
# Stair = getStairCase(12)
# @profview res = SW.getAllNeighborStates(Stair)
##

function solveED(state,μ,args...;kwargs...)
    @time res = SW.getAllNeighborStates(state)
    @time H = SW.H(res.AllStates,res.Neighbors,μ)
    @time sol = SW.SolveH(H,args...;kwargs...)
    return (;res,sol,H)
end

function getObservables(res,sol)
    # state = res.AllStates[1].parent

    AllStates = fetch.([Threads.@spawn SW.spinConfig(state) for state in res.AllStates ])

    m = SW.getMagnetization(AllStates,sol)

    mavg = sum(abs,m)/length(m)
    mavgBulk = sum(abs,m[3:end-2,3:end-2])/length(m[3:end-2,3:end-2])
    magnetization = (;m,mavg,mavgBulk)


    structureFac = SW.getStructureFac(AllStates,sol)

    return (;magnetization...,structureFac...)

end
##
function plotMagnetization!(ax,observables)
    heatmap!(ax,observables.m,colormap = :grays)
end

function plotStructureFac!(ax,observables)
    (;k,Sq_k) = observables
    heatmap!(ax,k,k,real(Sq_k))
end

function plotOverview(Observables;title=L"")
    with_theme(theme_SimpleTicks()) do  
        fig = Figure(;size = 0.8 .*(450,600))

        axmag = Axis(fig[1,1];title,xlabel = L"x",ylabel = L"y",aspect = 1)
        axSq = Axis(fig[2,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,xticks = PiTicks(0:0.5pi:2pi),yticks = PiTicks(0:0.5pi:2pi))
        
        hm = plotMagnetization!(axmag,Observables)
        Colorbar(fig[1,2],hm,label = L"\langle S^z \rangle")
        hm = plotStructureFac!(axSq,Observables)
        Colorbar(fig[2,2],hm,label = L"\mathcal{S}^{zz}(\mathbf{q})")
        fig
    end

end

plotOverview(res,sol;kwargs...) = plotOverview(getObservables(res,sol);kwargs...)

##

μ = 0
Sol5x5 = solveED(SW.periodicState5x5(14),μ)
##
# @profview plotOverview(Sol5x5.res,Sol5x5.sol,title = L"$5×5$ state, $μ = %$μ$")
@time plotOverview(Sol5x5.res,Sol5x5.sol,title = L"$5×5$ state, $μ = %$μ$")
##

Sol6x6 = solveED(SW.periodicState6x6(14),μ)
@time plotOverview(Sol6x6.res,Sol6x6.sol,title = L"$6×6$ state, $μ = %$μ$")
##
SolStair = solveED(SW.getStairCase(14),μ)
##
@time plotOverview(SolStair.res,SolStair.sol,title = L"stair state, $μ = %$μ$")
##

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

StairEnergy(L) = getEnergy(getStairCase(L))
# Energy5x5(L) = -L^2/10
# Energy5x5(L) = getEnergy(SW.periodicState5x5(L))
Energy5x5(L) = -length(SW.getApplicablePlaquettes(SW.periodicState5x5(L)))
Energy6x6(L) = getEnergy(SW.periodicState6x6(L))
Energy6x6_3(L) = getEnergy(SW.periodicState6x6_3(L))
Ls_small = 5:1:13
Ls_medium = 5:1:13
Ls_large = 5:1:13
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
    hlines!(ax,1/10,label = L"E_0 =-L^2/10",linestyle = :dash, color = :grey)
    hlines!(ax,1/16,label = L"E_0 =-L^2/16",linestyle = :dashdot, color = :grey)
    axislegend(ax,position = :rb,nbanks = 2)
    save("exactFig/EnergyScaling.png",fig)
    fig

end
##_______________________Periodic Boundaries_______________________

Stair = getPeriodic(getStairCase(8))
SW.plotApplPlaquettes(Stair)
##

SW.flipPlaquette!(Stair,6,1)
SW.plotApplPlaquettes(Stair)

##

function plotSpectrum!(ax,axgap,L;kwargs...)
    State = getStairCase(L)
    # Stair = SW.periodicState5x5(L) 
    eval = solveED(State,0,nev = 15).sol.values
    Δ = eval[2] - eval[1]
    scatter!(ax,1:length(eval),eval;kwargs...)
    scatter!(axgap,L,Δ;kwargs...)
end

let 
    fig = Figure(size = 400 .* (1,1.5))
    axgap = Axis(fig[1,1],xlabel = L"L",ylabel = L"Δ")
    ax = Axis(fig[2,1],xlabel = L"n",ylabel = L"E_n")
    plotSpectrum!(ax,axgap,8,markersize = 10,label = L"L=8")
    plotSpectrum!(ax,axgap,10,markersize = 10;color = :blue,label = L"L=10")
    plotSpectrum!(ax,axgap,12,markersize = 10;color = :red,label = L"L=12")
    plotSpectrum!(ax,axgap,13,markersize = 10;color = :darkred,label = L"L=13")
    plotSpectrum!(ax,axgap,14,markersize = 10;color = :black,label = L"L=14")
    rowsize!(fig.layout,1,Relative(0.4))
    axislegend(ax,position = :lt,nbanks = 2)
    fig
end
##

StairEnergy(L) = getEnergy( getPeriodic(getStairCase(L)))
Energy5x5(L) = getEnergy( getPeriodic(SW.periodicState5x5(L)))
# Energy5x5(L) = -length(SW.getApplicablePlaquettes(getPeriodic(SW.periodicState5x5(L))))
Energy6x6(L) = getEnergy( getPeriodic(SW.periodicState6x6(L)))
Energy6x6_3(L) = getEnergy( getPeriodic(SW.periodicState6x6_3(L)))
Ls_4 = 4:4:12
Ls_5 = 5:5:15
Ls_6 = 6:6:12
StairEs = ([StairEnergy(L) for L in Ls_4])
Ens_5x5 = ([Energy5x5(L) for L in Ls_5])
Ens_6x6 = ([Energy6x6(L) for L in Ls_6])
Ens_6x6_3 = ([Energy6x6_3(L) for L in Ls_6])

##
with_theme(theme_SimpleTicks()) do
    fig = Figure(size = 0.8 .*(600,400))
    ax = Axis(fig[1,1],xlabel = L"L",ylabel = L"-E_0(L)/L^2",
    xscale = log10,yscale = log10
    )
    scatterlines!(ax,Ls_4,-StairEs ./Ls_4.^2 ,label = L"Staircase$$",color = :red)
    scatterlines!(ax,Ls_5,-Ens_5x5 ./Ls_5.^2 ,label = L"5×5",color = :blue)
    scatterlines!(ax,Ls_6,-Ens_6x6 ./Ls_6.^2 ,label = L"6×6",color = :green)
    scatterlines!(ax,Ls_6,-Ens_6x6_3 ./Ls_6.^2 ,label = L"⋱_{6×6}", color = :black)
    LRange = LinRange(4,15,100)
    lines!(ax,LRange,1/10 .*LRange.^0 ,label = L"E_0 =-L^2/10",linestyle = :dash, color = :grey)
    lines!(ax,LRange,1/16 .*LRange.^0 ,label = L"E_0 =-L^2/16",linestyle = :dashdot, color = :grey)
    axislegend(ax,position = :lt,nbanks = 2)
    ylims!(ax,0.05,0.25)
    save("exactFig/EnergyScaling_periodic.png",fig)
    fig

end
##
Sol5x5 = solveED(getPeriodic(SW.periodicState5x5(15)),μ)
##
@time plotOverview(Sol5x5.res,Sol5x5.sol,title = L"$5×5$ state, $μ = %$μ$")

##
# Hardcore bosons
function isAllowedState(state)
    for i in axes(state,1),j in axes(state,2)
        if state[i,j]
            NearestNeighborOccupied = state[i-1,j] || state[i+1,j] || state[i,j-1] || state[i,j+1]

            if NearestNeighborOccupied
                return false
            end

            NextNearestNeighborOccupied = state[i-1,j-1] || state[i+1,j-1] || state[i-1,j+1] || state[i+1,j+1]

            if NextNearestNeighborOccupied
                return false
            end
        end
    end
    return true
    
end

function getAllStates(L)
    states = 0:1
    combinations = Iterators.product((states for i in 1:L, j in 1:L)...)

    allowedStates = BitMatrix[]

    state = SW.PeriodicMatrix(BitMatrix(zeros(Bool,L,L)),L,L,0)
    for c in combinations

        for i in eachindex(state.UC)
            state.UC[i] = c[i]
        end

        if isAllowedState(state)
            push!(allowedStates,copy(state))
        end
    end
    return allowedStates

end

@time getAllStates(4)
