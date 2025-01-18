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


with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"\mu",ylabel = L"E/L^2")
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
##
Sqs = stack([SW.normalized_Sq(h5read(file,"Sq_numerator"),h5read(file,"obs_denominator")) for file in readdir("../../Data/obsS12_staircase",join=true)])
Ens = stack([SW.normalized_En(h5read(file,"Energy"),h5read(file,"en_denominator")) for file in readdir("../../Data/obsS12_staircase",join=true)])


##
errlines(dropmean(Ens,dims=2),dropstd(Ens,dims=2))
##
with_theme(theme_SimpleTicks()) do 
    SqMat = dropmean(Sqs,dims=4)[:,:,5]
    SqErr = dropstd(Sqs,dims=4)[:,:,5]
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
    heatmap!(ax,qx,qy,Sq_q)
    fig
end
