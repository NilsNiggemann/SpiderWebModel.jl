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

## Spin 1


with_theme(theme_SimpleTicks()) do
    # Prepare configurations
    # 1. Single spin flip
    S1 = SW.stencilConfig(zeros(9,8), 1)
    S1[5,5] = 2

    # 2. Lineon
    S2 = SW.stencilConfig(zeros(9,8), 1)
    S2[5:2:end,5] .= 2

    # 3. Single fracton (flip one spin and plot fractons)
    S3 = SW.stencilConfig(zeros(9,8), 1)
    S3[5:2:end,1:2:5] .= 2

    # Create 1x4 figure
    fig = Figure(size = (4*160, 160), backgroundcolor = :transparent)
    configs = [S1, S2, S3]

    for (i, S) in enumerate(configs)
        ax = Axis(fig[1, i];
            SW.getConfigAxis(S)...,
            xticks = SimpleTicks(1:10), yticks = 1:10,
            xticksvisible = false, yticksvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false,
            backgroundcolor = :transparent,
        )
        SW.plotFractons!(ax, S, markersize = 19)
    end

    elem_1 = [MarkerElement(color = :black, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    elem_2 = [MarkerElement(color = :grey, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    elem_3 = [MarkerElement(color = :white, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]

    Legend(fig[1, 3, Right()],
    [elem_1,elem_2,elem_3, ],
    [L"|↓⟩",L"|0⟩",L"|↑⟩"],
    patchsize = (35, 35), rowgap = 10,nbanks=1,backgroundcolor = (:black,0.05),colgap=1,padding = (2,2,-2,-2),patchlabelgap=-2)

    save("Fracton_overview.svg", fig)
    fig
end
