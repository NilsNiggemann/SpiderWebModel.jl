import Pkg

import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
meanstd(x) = (mean(x),std(x))
include("4x4Orders_base.jl")
##

##
a = generatePeriodic(4,1)
a_st = sort!(collect(filterConfs(a,1)),by=x->sum(abs,x))


a_configs = let
    L = 4
    confs = [makeConf(UC,L,1) for UC in a_st]

    # filter!(SW.fulFillsConstraint, confs)
    filter!(x->length(SW.getApplicablePlaquettes(x)) > 0,confs)
    @assert all(SW.fulFillsConstraint,confs)
    confs
end
##
reducedConfigs = getMaxFlipConfs(a_configs;Nwalkers = 28,NSteps = 1) #first reduction
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 10)
##
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28,NSteps = 100)
##
reducedConfigs = makeConf.(filterConfs(parent.(parent.(reducedConfigs)),1),4,1)
reducedConfigs = getMaxFlipConfs(reducedConfigs;Nwalkers = 28*1,NSteps = 1000)

##
SW.h5write("Data/reducedConfigs.h5","reducedConfigs",stack(reducedConfigs))

##
L = 12
reducedConfigs = makeConf.(collect.(eachslice(SW.h5read("Data/reducedConfigs.h5","reducedConfigs"),dims=3)),4,1)
AllSectors = makeConf.(reducedConfigs,L,1)
##

function plotConfs!(fig,Confs; transpose = false)
    numGS = length(Confs)

    ijInds = [(i,j) for i in 1:round(Int,sqrt(numGS),RoundUp),j in 1:round(Int,sqrt(numGS))]
    if transpose
        ijInds = [(j,i) for i in 1:round(Int,sqrt(numGS),RoundUp),j in 1:round(Int,sqrt(numGS))]
    end
    
    with_theme(theme_SimpleTicks()) do

        axs = [Axis(fig[i,j]; SW.getConfigAxis(Confs[1])...,xlabelvisible=false,ylabelvisible=false,xticklabelsvisible=false,yticklabelsvisible=false,aspect=1,xticks = 1:4,yticks = 1:4) for (i,j) in ijInds[1:length(Confs)]]

        for ind  in eachindex(Confs)
            ax = axs[ind]
            # newConf = similar(Confs[i],4,4)
            # newConf[1:4,1:4] .= @view Confs[i][1:4,1:4]
            SW.plotApplPlaquettes!(ax,Confs[ind])

            i,j = ijInds[ind]
            Label(fig[i,j, TopLeft()],L"%$ind$$",fontsize = 12,padding = (-15,0,-10,0),color = :black)
            # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-20,0))

            # text!(ax,Point(1,4),text = "$i",color = :lime,fontsize = 15,align = (:center,:center))

        end
        # unique_i = unique(i for (i,j) in ijInds)
        # unique_j = unique(j for (i,j) in ijInds)
        # for i in unique_i[1:end-1]
        #     rowgap!(fig.layout,i,0.)
        # end
        # for j in unique_j[1:end-1]
        #     colgap!(fig.layout,j,10)
        # end
        return fig
    end 

end
function plotConfs(Confs;kwargs...)
    fig = Figure(fontsize = 22,size = length(Confs) .*(12,12))
    plotConfs!(fig,Confs;kwargs...)
end

plotConfs(reducedConfigs)
##
μ = 0.95

CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(AllSectors[1]),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))

##
@time res = findEnergies(AllSectors,CT,ψG;Nwalkers = 20,NSteps = 800)
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"# Config $$",ylabel = L"E_0/N_\text{sites}")
    Nsites = 4*4
    errorbars!(eachindex(res.en),res.en ./Nsites,res.Δen ./Nsites,whiskerwidth=4,color = :black)
    scatter!(res.en ./Nsites,marker = '×',color = :black)
    fig
end
##
Allres = empty!([res])
mus_sectors = LinRange(-0.1,0.99,20)
for mu in mus_sectors
    CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(AllSectors[1]),SW.Hxx_RK(mu))
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-mu))
    res = findEnergies(AllSectors,CT,ψG;Nwalkers = 28*6,NSteps = 3000)
    push!(Allres,res)
end

##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 26,size = 2 .* (700,400))
    E_scal = 0.2
    ens = getproperty.(Allres,:en)
    Δens = getproperty.(Allres,:Δen)
    mus = mus_sectors
    Elin = E_scal .* (mus .-1)

    ax = Axis(fig[1,1],xlabel = L"\mu",ylabel = L"E_0/(N_\text{sites}(1 - \mu ))")
    Nsites = 12*12

 
    minsectors = argmin.(ens)

    text_locations = Point2f[]

    linewidths = LinRange(1.2,3.0, length(reducedConfigs))

    function getLS(width)
        if width < 2
            return :solid
        else
            return :dash
        end
    end
    # Generate all combinations of linestyles and linewidths

    for i_sector in eachindex(reducedConfigs)
        linewidth = linewidths[i_sector]
        linestyle = getLS(linewidth)
        
        en = getindex.(ens,i_sector) ./ Nsites ./(1 .-mus)
        Δen = getindex.(Δens,i_sector)./(1 .-mus) ./ Nsites
        l = scatterlines!(ax,mus,en ,marker = '×';linestyle,linewidth,markersize = 15)
        errorbars!(ax,mus,en ,Δen,whiskerwidth=4)
        # text!(Point(-0.15,en[1]),text = "%$i_sector")

        # pos = round(Int,length(en)/length(AllSectors) * i_sector)
        # pos = max(1,pos)
        # pos = min(length(en),pos)

        # text!(ax,Point(mus[pos],en[pos]-0.01),text = L"%$i_sector", color = l.color[],fontsize = 10)
        # ppoint = isodd(i_sector) ? Point(mus[1]-0.05,en[begin]) : Point(mus[end]+0.05,en[end])
        ppoint = Point2f(mus[1]-0.03,en[begin])
        ppointRd = Point2f(ppoint[1],round(ppoint[2],digits = 2))
        
        while ppointRd in text_locations
            ppoint = Point2f(ppoint[1]+0.02,ppoint[2])
            ppointRd = Point2f(ppoint[1],round(ppoint[2],digits = 2))
        end
        push!(text_locations,ppointRd)
        
        text!(ax,ppoint,text = L"%$i_sector", color = l.color[],fontsize = 12,align = (:top,:top))
    end

    for i in eachindex(mus)[3:end]
        if minsectors[i-1] != minsectors[i]
            vlines!(ax,mus[i-1],color = :black,linestyle = :dash,linewidth = 0.8)
        end
    end
    colsize!(fig.layout,1,Relative(0.6))
    newfig = plotConfs!(fig[:,2:3],reducedConfigs,transpose=true)
    
    confaxs = size(newfig.layout)
    for i in 1:confaxs[1] -1
        rowgap!(newfig.layout,i,0.)
    end
    for i in 2:confaxs[2]-1
        colgap!(newfig.layout,i,0.)
    end
    # fig_confs = GridLayout() = newfig
    # fig.layout[1,2] = fig_confs
    # ax2 = Axis(fig[1,1],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    # ax3 = Axis(fig[1,2],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    # ax4 = Axis(fig[1,3],aspect=1,xticks = 1:2:4,yticks = 1:2:4)
    # rowsize!(fig.layout,1,Relative(0.3))
    # SW.plotApplPlaquettes!(ax2,reducedConfigs[1])
    # SW.plotApplPlaquettes!(ax3,reducedConfigs[5])
    # SW.plotApplPlaquettes!(ax4,reducedConfigs[6])
    fig
end
##
# perm = sortperm(res.en)
# # SW.plotApplPlaquettes(a_configs[perm][3])

# emin = minimum(res.en)
# numGS = findfirst(>(emin+3e-1),res.en[perm])

##

##
μ = -0.1

CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(a_configs[1]),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))


@time results1 = [SW.startManyWalkerGFMC(upscale(reducedConfigs[end-5],8),CT,28*6,5000,ψG) for _ in 1:10]
##
@time results2 = [SW.startManyWalkerGFMC(upscale(reducedConfigs[end-6],8),CT,28*6,5000,ψG) for _ in 1:10]
##
plotEnergies(results1,CT;normalize=true,dense=true,τ = 10,color = :black)
plotEnergies!(results2,CT;normalize=true,dense=true,τ = 10,color = :red)
current_figure()

##
S = upscale(reducedConfigs[1],28)

μ = 0.4
CT = SW.ContinuousTimeMethod(0.1,1,-0.266length(S),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
SW.plotApplPlaquettes(S)
##
results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,28*30,10000,ψG) for _ in 1:1])
##
plotEnergies(results,CT;normalize=true,dense=true,τ = 20,color = :black)
##
Sqs = SW.getSqsGFMC(results, round(Int,10 ÷ CT.τ);nBra = 1)
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    ax2 = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    SqMat = mean(Sqs)
    SqErr = std(Sqs)
    fittingCoefs = optimizeCoeffs(SqMat)
    Sq = SW.getSqCont(SqMat)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    heatmap!(ax,qx,qy, Sq.(qs)
    )
    SQFT(x) = SqFieldTheory(SVector(x),fittingCoefs)
    heatmap!(ax2,qx,qy,SQFT.(qs))
    fig
end