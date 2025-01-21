import Pkg
cd(@__DIR__)
Pkg.activate("../../")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using SpiderWebModel.HDF5
include("../plottingUtils.jl")

muRange = h5read("../../Data/energy_mu_S12.h5","muRange")
ens1 = h5read("../../Data/energy_mu_S12.h5","energy1")
ens2 = h5read("../../Data/energy_mu_S12.h5","energy2")
ens3 = h5read("../../Data/energy_mu_S12.h5","energy3")
ens4 = h5read("../../Data/energy_mu_S12.h5","energy4")

ens = [ens1,ens2,ens3,ens4]
E_trivial(mu) = -0.1*(1-mu)
function getSectorConfig(L,i)
    S = SW.stencilConfig(0.5*ones(L,L),1/2;boundaryCondition = :periodic)

    S.= SW.getSelectedS12PeriodicState(L,i)
end

spinecolors(color) = (;topspinecolor = color,bottomspinecolor = color,leftspinecolor = color,rightspinecolor = color)

function plotSectorEnergies!(fig,ax,ens,configs;inset_scale = 1,top_inset_pos=idx->(0.05 +0.18(idx-1),-0.01),bot_inset_pos=(0.9,-0.09),axkwargs...)
    axkwargs = SW.getConfigAxis(getSectorConfig(8,1))
    colors = [:blue,:green,:red,:purple]

    inax = [
        insetAtPoint(fig,ax,top_inset_pos(idx),inset_scale.*(36,36);
        spinecolors(colors[i])...,
        spinewidth = 4,
        xticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        axkwargs...
        ) for (idx,i) in enumerate(configs)
    ]

    inaxTriv = insetAtPoint(fig,ax,bot_inset_pos,inset_scale.*(50,50);
        spinecolors(:grey)...,
        spinewidth = 4,
        xticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        title = L"Trivial $$",
        axkwargs...
    )

    SW.plotApplPlaquettes!(inaxTriv,getSectorConfig(20,5),markersize = 8)

    etriv = E_trivial.(muRange)# ./ (1 .-muRange)

    lines!(ax,muRange,etriv,color = :grey,linewidth = 3)
    # ylims!(ax,-0.12,0.03)
    Ls = [20,24,20,18]
    for (idx,i) in enumerate(configs)
        SW.plotApplPlaquettes!(inax[idx],getSectorConfig(Ls[i],i),markersize = 8)
    end

    for (i,en) in enumerate(ens)
        L = Ls[i]
        color = colors[i]
        e = dropmean(en,dims=2) ./ L^2 # ./ (1 .-muRange)
        e_err = dropstd(en,dims=2) ./ L^2 # ./ (1 .-muRange)
        scatterlines!(ax,muRange,e;color,marker = '×',markersize = 15)
        errorbars!(ax,muRange,e,e_err;color)
        # errlines!(ax,muRange,e,e_err;color)
    end
    fig
end

with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"\mu",ylabel = L"E/L^2")


    plotSectorEnergies!(fig,ax,ens,[1,2,3,4])
end
##
Sqs_mu0 = stack([h5read(file,"StructureFactor") for file in readdir("../../Data/obsS12_staircase",join=true)])
Ens = stack([h5read(file,"Energy") for file in readdir("../../Data/obsS12_staircase",join=true)])


##
errlines(dropmean(Ens,dims=2),dropstd(Ens,dims=2))
##
with_theme(theme_SimpleTicks()) do 
    SqMat = dropmean(Sqs_mu0,dims=4)[:,:,end]
    SqErr = dropstd(Sqs_mu0,dims=4)[:,:,end]
    # SqMat = Sqs[:,:,5]
    # SqErr = Sqs[:,:,5]
    fig = Figure(size = 120 .* (4,4),fontsize = 22)

    xticks = yticks = PiTicks([0,pi])
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = true,yticklabelsvisible = true)

    # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
    Sq = SW.getSqCont(SqMat)
    Sqerr = SW.getSqCont(SqErr)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    hm = heatmap!(ax,qx,qy,Sq_q)
    Colorbar(fig[1,2],hm,label = L"\langle S^z(\mathbf{q})^*  S^z(\mathbf{q})\rangle")

    fig
end

##
ClassicalConfigs = let 
    files = readdir("../../ConfsRaw/GurobiConfs",join=true)

    Confs = Vector{Matrix{Bool}}()

    for file in files
        h5open(file,"r") do f
            for k in keys(f["Confs"])
                confs = read(f["Confs/$k"])
                for conf in eachslice(confs,dims=3)
                    push!(Confs, collect(conf))
                end
            end
        end
    end
    Spinconfs = SW.floatSpinConfig.(Confs,0.5)
    @assert all(SW.fulFillsConstraint,Spinconfs)
    Spinconfs
    
end
SqClassical = SW.getEqualWeightStructureFac(ClassicalConfigs)
Sq_RK = h5read("../../Data/StaircaseS12_L40_RK.h5","Sqs") ./4
##
function sumRuleCheck(SqMat)
    finiteEls = filter!(isfinite,SqMat[1:end-1,1:end-1][:])
    return sum(finiteEls)/length(finiteEls)
end
##
function makeSqPlot(Sq_Classical,Sq_mu0,Sq_RK)
        


    fig = Figure(size = 1.8 .*(520,300),fontsize=20)
    ticks = PiTicks([0,pi])

    L40RescalingFactor = 3

    Left_Half = GridLayout()
    Right_Half = GridLayout()
    SpinConf_Fig = GridLayout()
    StrucFac_Top = GridLayout()
    StrucFac_Bot = GridLayout()
    EN_Plot = GridLayout()

    fig.layout[1:2,1] = Left_Half
    fig.layout[1:2,2] = Right_Half
    
    Left_Half[1, 1] = StrucFac_Top
    Left_Half[2, 1] = StrucFac_Bot
    Right_Half[1, 1] = SpinConf_Fig
    Right_Half[2, 1] = EN_Plot

    exemplaryConfig = SW.SpinConfig(ClassicalConfigs[1][10:22,10:22],0.5)
    stairCase_state = SW.SpinConfig(SW.getStairCase(12),0.5)

    axkwargs = (;xminorticksvisible=true,yminorticksvisible=true,xtickwidth = 1.4,ytickwidth = 1.4,xminortickwidth = 1.2,yminortickwidth = 1.2,xticksize=5,yticksize=5,xminorticksize=3,yminorticksize=3,xminortickalign =1,yminortickalign =1,xtickalign =1,xticksmirrored =true,yticksmirrored=true,ytickalign=1,
    xlabelpadding = -5,
    xtickcolor = :white,
    xminortickcolor = :white,
    ytickcolor = :white,
    yminortickcolor = :white,
    
    )


    SpinConf_Fig[1, 1] = ax_generic_conf = Axis(fig; SW.getConfigAxis(exemplaryConfig)...,xticklabelsvisible = false,yticklabelsvisible = false
    )

    SW.plotApplPlaquettes!(ax_generic_conf,exemplaryConfig,markersize = 8)

    SpinConf_Fig[1, 2] = ax_staircase_conf = Axis(fig; SW.getConfigAxis(stairCase_state)...,xticklabelsvisible = false,yticklabelsvisible = false,spinecolors(:blue)...,spinewidth = 5
    )
    SW.plotApplPlaquettes!(ax_staircase_conf,stairCase_state,markersize = 8)

    StrucFac_Top[1, 1] = ax = Axis(fig, aspect = 1,xticks = ticks,yticks = ticks,xminorticksvisible = true ,xlabel = L"q_x",ylabel = L"q_y",yminorticksvisible = true,
    xlabelvisible=false,xticklabelsvisible=false;
    axkwargs...
    )

    StrucFac_Top[1, 2] = ax2 = Axis(fig, aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,
    xlabelvisible=false,xticklabelsvisible=false,xticks = ticks,yticks = ticks,
    # title = L"class. spin-$1/2$"
    ;
    axkwargs...
    )


    SqFunc = SW.getSqCont(real(Sq_Classical.Sq))

    kx = collect(trueMomenta(-0.5pi,1.5pi,size(Sq_Classical.Sq,1)-1))
    # kx = kx
    ky = kx
    # kxSpin1 = collect(trueMomenta(-0.5pi,1.5pi,size(SqSpin1,1)-1))
    kxSpin1 = collect(trueMomenta(-0.5pi,1.5pi,size(Sq_Classical.Sq,1)-1))
    # kxSpin1 = 2pi .* (0.5:1)
    kySpin1 = kxSpin1

    SqLN = [SqLargeN(kx,ky)/4 for kx in kxSpin1 , ky in kySpin1]
    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]
    @info "" sumRuleCheck(SqPlot) sumRuleCheck(SqLN)
    # cRange = extrema(SqPlot)
    cRange = (0,maximum(SqPlot))
    hm = heatmap!(ax,kx,ky,SqLN,colorrange = cRange)
    hm2 = heatmap!(ax2,kx,ky,SqPlot,colorrange = cRange)
    # hm = halfhalfheatmap!(ax,kx,ky,SqFunc,SqLargeN,x->-x+3pi,normalize = true)
    EN_Plot[1,1] = ax_sector_energies = with_theme(theme_SimpleTicks()) do 
        Axis(fig[1,2],xlabel = L"\mu",ylabel = L"E/L^2")
    end
    StrucFac_Bot[1, 1] = axmu0 = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
    axkwargs...,
    # xtickcolor = :white,ytickcolor = :white,xminortickcolor = :white,yminortickcolor = :white
    )
    StrucFac_Bot[1, 2] = axRK = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,ylabelvisible = false,yticklabelsvisible = false,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
    axkwargs...,
    # xtickcolor = :white,ytickcolor = :white,xminortickcolor = :white,yminortickcolor = :white
    )

    kx = ky = trueMomenta(-0.5pi,1.5pi,size(Sq_mu0,1)-1)

    SqFunc = SW.getSqCont(Sq_mu0 ./ L40RescalingFactor)
    SqMat = [real(SqFunc(x,y)) for x in kx, y in ky] 

    # SqRKFunc = SW.getSqCont(SqRK)
    # SqRKMat = [real(SqRKFunc(x,y)) for x in kx, y in ky]
    # colorrange = extrema(SqRKMat)

    
    kxGFMC = kyGFMC = trueMomenta(-0.5pi,1.5pi,size(Sq_RK,1)-1)
    
    Sq_RKFunc = SW.getSqCont(mean(eachslice(Sq_RK,dims=3)) ) 
    Sq_RK_Mat = [real(Sq_RKFunc(x,y)) for x in kxGFMC, y in kyGFMC]

    colorrange = (0,maximum(maximum.((Sq_RK_Mat,SqMat))))
    
    hm_RK = heatmap!(axRK, kx, ky, Sq_RK_Mat;colorrange,colormap = :turbo)
    hm_mu0 = heatmap!(axmu0, kx, ky, SqMat;colorrange,colormap = :turbo )



    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-20,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-20,0))
    # Label(fig[1,3, TopLeft()],L"c)$$",padding = (-30,0,-20,0))
    # Label(fig[2,1, TopLeft()],L"d)$$",padding = (-30,0,-20,0))
    # Label(fig[2,2, TopLeft()],L"e)$$",padding = (-30,0,-20,0))

    textpos = Point(-pi/2,3pi/2)


    Colorbar(StrucFac_Top[2,1:2], hm2,height = Relative(0.8),width = Relative(0.99),ticks = SimpleTicks([0,0.2,0.4,0.6]),
    vertical=false,
    # label = L"\mathcal{S}(q)"
    )
    Colorbar(StrucFac_Bot[0,1:2], hm_RK,height = Relative(0.8),width = Relative(0.99),ticks = SimpleTicks(),
    vertical = false, flipaxis = false,
    # label = L"\mathcal{S}(q)"
    )

    # Colorbar(StrucFac_Top[1,3], hm2,height = Relative(0.8),width = Relative(0.99),ticks = SimpleTicks([0,0.2,0.4,0.6]),
    # # label = L"\mathcal{S}(q)"
    # )
    # Colorbar(StrucFac_Bot[1,3], hm_RK,height = Relative(0.8),width = Relative(0.99),ticks = SimpleTicks(),
    # # label = L"\mathcal{S}(q)"
    # )

    
    # Label(StrucFac_Top[2,1:2, Top()],L"\textrm{Classical }\mathcal{S}(q)",padding = (0,0,20,0))
    # Label(StrucFac_Bot[0,1:2, Top()],L"\textrm{Quantum }\mathcal{S}(q)",padding = (0,0,-80,0))

    rowsize!(StrucFac_Top,2,Relative(0.06))
    rowsize!(StrucFac_Bot,0,Relative(0.06))

    colsize!(fig.layout,1,Relative(0.43))
    rowsize!(Right_Half,1,Relative(0.43))
    colgap!(fig.layout,1,5)

    rowgap!(Left_Half,1,0)

    colgap!(StrucFac_Top,1,0)
    rowgap!(StrucFac_Top,1,0)

    colgap!(StrucFac_Bot,1,0)
    rowgap!(StrucFac_Bot,1,0)

    # colgap!(StrucFac_Top,2,3)
    # colgap!(StrucFac_Bot,1,3)
    # colgap!(StrucFac_Bot,2,3)
    
    plotSectorEnergies!(fig,ax_sector_energies,ens,[2,3,4],inset_scale = 1,
    top_inset_pos=idx->(0.05 +0.195(idx-1),-0.015),bot_inset_pos=(0.9,-0.085)
    )
    
    
    sub_fig_labelsize = 22
    text!(ax, textpos ,text = L"a)",color = :black,align = (:left,:top),fontsize = sub_fig_labelsize)
    text!(ax2, textpos ,text = L"b)",color = :black,align = (:left,:top),fontsize = sub_fig_labelsize)

    Label(SpinConf_Fig[1, 1, TopLeft()], L"c)", padding = (-30, 0, -20, 0),fontsize = sub_fig_labelsize)
    Label(SpinConf_Fig[1, 2, TopLeft()], L"d)", padding = (-30, 0, -20, 0),fontsize = sub_fig_labelsize)
    Label(EN_Plot[1, 1, TopLeft()], L"g)", padding = (-30, 0, -20, 0),fontsize = sub_fig_labelsize)

    text!(axmu0, textpos ,text = L"e)",color = :white,align = (:left,:top),fontsize = sub_fig_labelsize)
    text!(axRK, textpos ,text = L"f)",color = :white,align = (:left,:top),fontsize = sub_fig_labelsize)
    if L40RescalingFactor != 1
        text!(axmu0, Point(pi,3pi/2) ,text = L"\times \frac{1}{%$(L40RescalingFactor)}",color = :white,align = (:left,:top),fontsize = sub_fig_labelsize)
    end

    # save("../../figs/Sq_comparison_2.png", fig,px_per_unit=4)
    # Colorbar(fig[1, 2], hm)
    fig
end
##
makeSqPlot(SqClassical,dropmean(Sqs_mu0,dims=4)[:,:,end],Sq_RK)
