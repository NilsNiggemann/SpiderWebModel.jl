import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
##

with_theme(theme_SimpleTicks()) do
    S = SW.stencilConfig(0.5 .* ones(9,8),1/2)
    S[5:2:end,1:2:5] .= -1

    # S[5:2:end,7] .*= -1
    fig = Figure(size = (20,20) .*size(S),backgroundcolor = :transparent)
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,xticks = SimpleTicks(1:10),yticks = 1:10,xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false,backgroundcolor = :transparent)
    SW.plotSpinConfig!(ax,S)
    SW.plotFractons!(ax,S,markersize = 19)
    save("fractons.svg",fig)
    fig
end
##
with_theme(theme_SimpleTicks()) do
    S = SW.stencilConfig(0.5 .* ones(9,8),1/2)

    S[5:2:end,5] .*= -1
    fig = Figure(size = (20,20) .*size(S),backgroundcolor = :transparent)
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,xticks = SimpleTicks(1:10),yticks = 1:10,xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false,backgroundcolor = :transparent)
    SW.plotSpinConfig!(ax,S)
    SW.plotFractons!(ax,S,markersize = 19)
    save("lineons.svg",fig)
    fig
end

##
with_theme(theme_SimpleTicks()) do
    # global i +=1
    SW.Random.seed!(34)
    S = SW.constructConfigPath(6,6,SW.ALLGS_S12,)
    S = SW.SpinConfig(S[1:12,1:12],1/2)

    fig = Figure(size = (20,20) .*size(S))
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,xticks = SimpleTicks(1:10),yticks = 1:10,xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false)
    SW.plotSpinConfig!(ax,S)
    SW.plotApplPlaquettes!(ax,S)
    save("Application/Confs/example_2Plaquettes.svg",fig)
    fig
end


##
with_theme(theme_SimpleTicks()) do
    S = SW.SpinConfig(SW.PeriodicMatrix(SW.getStairCase(12)),1/2)
    fig = Figure(size = (20,20) .*size(S))
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,xticks = SimpleTicks(1:10),yticks = 1:10,xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false)
    SW.plotSpinConfig!(ax,S)
    SW.plotApplPlaquettes!(ax,S)
    save("Application/Confs/example_staircase.svg",fig)
    fig
end