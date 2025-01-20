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


function plotSectorEnergies!(fig,ax)
    axkwargs = SW.getConfigAxis(getSectorConfig(8,1))
    colors = [:blue,:green,:red,:purple]
    spincolors(color) = (topspinecolor = color,bottomspinecolor = color,leftspinecolor = color,rightspinecolor = color)

    inax = [
        insetAtPoint(fig,ax,(0.05 +0.18(i-1),-0.01),(36,36);
        spincolors(colors[i])...,
        spinewidth = 4,
        xticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        axkwargs...
        ) for i in eachindex(ens)
    ]

    inaxTriv = insetAtPoint(fig,ax,(0.9,-0.09),(50,50);
        spincolors(:grey)...,
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
    for (i,en) in enumerate(ens)
        L = Ls[i]
        SW.plotApplPlaquettes!(inax[i],getSectorConfig(Ls[i],i),markersize = 8)
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


    plotSectorEnergies!(fig,ax)
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
function makeSqPlot(Sq_Classical,Sq_mu0,Sq_RK)
        
    fig = Figure(size = 0.9 .*(420,340))
    ticks = PiTicks([0,pi])

    L40RescalingFactor = 3

    subgl_top = GridLayout()
    subgl_bot = GridLayout()

    axkwargs = (;xminorticksvisible=true,yminorticksvisible=true,xtickwidth = 1.4,ytickwidth = 1.4,xminortickwidth = 1.2,yminortickwidth = 1.2,xticksize=5,yticksize=5,xminorticksize=3,yminorticksize=3,xminortickalign =1,yminortickalign =1,xtickalign =1,xticksmirrored =true,yticksmirrored=true,ytickalign=1,
    xlabelpadding = -5,
    xtickcolor = :white,
    xminortickcolor = :white,
    ytickcolor = :white,
    yminortickcolor = :white,
    
    )

    subgl_top[1, 1] = ax = Axis(fig, aspect = 1,xticks = ticks,yticks = ticks,xminorticksvisible = true ,xlabel = L"q_x",ylabel = L"q_y",yminorticksvisible = true,
    xlabelvisible=false,xticklabelsvisible=false;
    axkwargs...
    )
    subgl_top[1, 2] = ax2 = Axis(fig, aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,
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
    @info "" sum(SqLN)/length(SqLN) sum(SqPlot)/length(SqPlot)
    # cRange = extrema(SqPlot)
    cRange = (0,maximum(SqPlot))
    hm = heatmap!(ax,kx,ky,SqLN,colorrange = cRange)
    hm2 = heatmap!(ax2,kx,ky,SqPlot,colorrange = cRange)
    # hm = halfhalfheatmap!(ax,kx,ky,SqFunc,SqLargeN,x->-x+3pi,normalize = true)

    subgl_bot[1, 1] = axmu0 = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
    axkwargs...,
    # xtickcolor = :white,ytickcolor = :white,xminortickcolor = :white,yminortickcolor = :white
    )
    subgl_bot[1, 2] = axRK = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,ylabelvisible = false,yticklabelsvisible = false,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
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
    
    hmED2 = heatmap!(axRK, kx, ky, Sq_RK_Mat;colorrange)
    hmED1 = heatmap!(axmu0, kx, ky, SqMat;colorrange )

    fig.layout[1, 1] = subgl_top
    fig.layout[2, 1] = subgl_bot
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-20,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-20,0))
    # Label(fig[1,3, TopLeft()],L"c)$$",padding = (-30,0,-20,0))
    # Label(fig[2,1, TopLeft()],L"d)$$",padding = (-30,0,-20,0))
    # Label(fig[2,2, TopLeft()],L"e)$$",padding = (-30,0,-20,0))

    textpos = Point(-pi/2,3pi/2)


    Colorbar(subgl_top[1,3], hm2,height = Relative(0.95),width = Relative(0.8),label = L"\mathcal{S}(q)",ticks = SimpleTicks([0,0.2,0.4,0.6]))
    Colorbar(subgl_bot[1,3], hmED2,height = Relative(0.95),width = Relative(0.8),label = L"\mathcal{S}(q)",ticks = SimpleTicks())


    text!(ax, textpos ,text = L"a)",color = :black,align = (:left,:top),fontsize = 18)
    text!(ax2, textpos ,text = L"b)",color = :black,align = (:left,:top),fontsize = 18)

    text!(axmu0, textpos ,text = L"c)",color = :white,align = (:left,:top),fontsize = 18)
    text!(axRK, textpos ,text = L"d)",color = :white,align = (:left,:top),fontsize = 18)
    text!(axmu0, Point(pi,3pi/2) ,text = L"\times \frac{1}{%$(L40RescalingFactor)}",color = :white,align = (:left,:top),fontsize = 14)
    
    colsize!(subgl_top,3,Relative(0.05))
    colsize!(subgl_bot,3,Relative(0.05))
    rowgap!(fig.layout,1,3)
    colgap!(subgl_top,1,3)
    colgap!(subgl_top,2,3)
    colgap!(subgl_bot,1,3)
    colgap!(subgl_bot,2,3)

    # ax_sector_energies = fig[1,4]
    # plotSectorEnergies!(fig,ax_sector_energies)
    # save("../../figs/Sq_comparison_2.png", fig,px_per_unit=4)
    # Colorbar(fig[1, 2], hm)
    fig
end
##
makeSqPlot(SqClassical,dropmean(Sqs_mu0,dims=4)[:,:,end],Sq_RK)
