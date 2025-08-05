import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
##
with_theme(theme_SimpleTicks()) do
    # Prepare configurations
    # 1. Single spin flip

    M⎸ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongCol!(M⎸,3,x->x-2)
    
    M_ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongRow!(M_,3,x->x-2)

    M╱ = SW.stencilConfig(zeros(8,8),1)

    # diag = SW.getDiagonal(S,3,1,true)
    # transformSpins!(diag,1)

    SW.transformSpinsAlongDiagonal!(M╱,3,1,x->x-2)

    M╲ = SW.stencilConfig(zeros(8,8),1)
    SW.transformSpinsAlongDiagonal!(M╲,5,-1,x->x-2)

    Confs = [M⎸,M_,M╱,M╲]

    # Create 1x4 figure
    # fig = Figure(size = (4*160, 160), backgroundcolor = :transparent)
    fig = Figure(backgroundcolor = :transparent)

    M = SW.stencilConfig(zeros(8,8),1)
    col = CartesianIndices(M)[2,1:2:end]

    M_ = SW.stencilConfig(zeros(8,8),1)
    row = CartesianIndices(M)[1:2:end,2]

    diaginds = SW.getDiagInds(M,2,1,true)

    antidiaginds = SW.getDiagInds(M,6,-1,true)

    changesites = [col,row,CartesianIndices(M)[diaginds],CartesianIndices(M)[antidiaginds]]
    # labels = ["M⎸","M_","M╱","M╲"]
    # $M_{\bm \vert}$, $M_{\bm -}$, $M_{\bm \diagup}$, and $M_{\bm \diagdown}$ 
    labels = [
        L"M_{\mathbf{|}} = \sum_{i ∈ □} {S^z_{i}}", 
        L"M_{\mathbf{-}} = \sum_{i ∈ □} {S^z_{i}}",
        L"M_{\mathbf{\diagup}} = \sum_{i ∈ □} {S^z_{i}}", 
        L"M_{\mathbf{\diagdown}} = \sum_{i ∈ □} {S^z_{i}}"
        ]
        
    colors = [:red,:blue,:lime,:orange]
    ax = Axis(fig[1,1];
    SW.getConfigAxis(M)...,
    backgroundcolor=:transparent,
    xticks = SimpleTicks([2,4,6,8]),
    yticks = SimpleTicks([2,4,6,8]),
    # xticklabelsvisible=false,
    # yticklabelsvisible=false,
    )
    ax2 = Axis(fig[1,2];
    SW.getConfigAxis(M)...,
    backgroundcolor=:transparent,
    xticks = SimpleTicks([2,4,6,8]),
    yticks = SimpleTicks([2,4,6,8]),
    # xticklabelsvisible=false,
    yticklabelsvisible=false,
    )
    SW.plotSpinConfig!(ax,M)
    SW.plotSpinConfig!(ax2,M)

    bound_points = 0.5*[
        0.999999*Point2f(-1,-1),
        0.999999*Point2f(1,-1),
        0.999999*Point2f(1,1),
        0.999999*Point2f(-1,1),
        0.999999*Point2f(-1,-1),
    ]


    col_alpha = 0.2
    function outline_site!(ax,site,args...;kwargs...)
        Ps = [Point(Tuple(site))-bp for bp in bound_points]
        lines!(ax,Ps,args...;kwargs...)
    end

    for i in 1:2
        for s in changesites[i]
            outline_site!(ax,s;color=colors[i],linewidth = 2)
            SW.plotSiteHighlight!(ax,Tuple(s);color = (colors[i],col_alpha))
        end
    end

    for i in 3:4
        for s in changesites[i]
            outline_site!(ax2,s;color=colors[i],linewidth = 2)
            SW.plotSiteHighlight!(ax2,Tuple(s);color = (colors[i],col_alpha))
        end
    end
    # elem_1 = [MarkerElement(color = :black, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    # elem_2 = [MarkerElement(color = :grey, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    # elem_3 = [MarkerElement(color = :white, marker = '◼', markersize = 30,strokecolor = :black,strokewidth=2)]
    elems = [[MarkerElement(color = (color,col_alpha), marker = '◼', markersize = 30,strokecolor = color,strokewidth=2)] for color in colors]

    Legend(fig[1, 2, Right()],
    elems,
    labels,
    patchsize = (35, 35), rowgap = 10,nbanks=1,backgroundcolor = (:black,0.05),colgap=1,padding = (2,2,-2,-2),patchlabelgap=-2)
    save("../figs/PaperFigs/ConservationLaws.svg",fig)

    fig
end
