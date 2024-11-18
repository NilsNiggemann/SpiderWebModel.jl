using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
include("../plottingUtils.jl")
function getRes(folder)
    files = let
        filesunsrt = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
        mus = [h5read(file,"mu") for file in filesunsrt]
        filesunsrt[sortperm(mus)]
    end
        #[8:end]
    # files = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=30/") for file in files]#[1:2:end]
    # filter!(!contains("mu=-0.18"),files)
    # filter!(!contains("mu=-0.2"),files)
    ##
    energies = stack([h5read(file,"energies") for file in files])
    mus = [h5read(file,"mu") for file in files]
    Sqs = stack([h5read(file,"SqsGFMC/100") for file in files])
    taus = [h5read(file,"tau") for file in files]
    return (;energies,mus,Sqs,taus,files)
end
##
outfilesFolder = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/"

res = getRes(outfilesFolder)
function filterRes(res,key)

    inds = findall(contains(key),res.files)

    return @views (;energies = res.energies[:,:,inds],mus = res.mus[inds],Sqs = res.Sqs[:,:,:,inds],taus = res.taus[inds],files = res.files[inds])
end
##
function getParentState(SECTOR_NAME,L)
    
    S0 = SW.stencilConfig(
        zeros(L,L),1,
        boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    if SECTOR_NAME == "S0"
    elseif SECTOR_NAME == "DenseLoops"
        S0 .= SW.periodicStateDenseLoops(size(S0,1))
    elseif SECTOR_NAME == "Diag"
        S0 .= SW.periodicStateDiag(size(S0,1))
    end
    S0
end

##
with_theme(theme_SimpleTicks()) do
    en_slice = res.energies[50,:,:]
    mus = res.mus
    fig = Figure(size = (450, 400))
    ax = Axis(fig[1, 1:3], xlabel = L"\mu", ylabel = L"E_0 / N_\textrm{sites}",xticks = SimpleTicks(-0.25:0.25:1.0),yticks = SimpleTicks(-0.3:0.1:0))

    sectors = [
    "S0",
    "DenseLoops",
        "Diag",
    ]
    labels = [
        L"S^z=0",
        L"Dense Loops $$",
        L"Diag $$",
    ]
    
    colors = [:black,:blue,:red]
    for (label,s,color) in zip(labels,sectors,colors)
        resfilt = filterRes(res,s)
        Nsites = prod(size(resfilt.Sqs)[1:2])
        emean = dropdims(mean(resfilt.energies[50,:,:],dims=1),dims=1) ./ Nsites
        # return emean, resfilt.mus
        enstd = dropdims(std(resfilt.energies[50,:,:],dims=1),dims=1) ./ Nsites

        scatterlines!(ax,resfilt.mus,emean ./ (1 .-resfilt.mus);label,color,marker = 'x')
        errorbars!(ax,resfilt.mus,emean ./ (1 .-resfilt.mus),enstd ./ (1 .-resfilt.mus);color)
        # scatter!(ax,mus,en)
    end
    
    axislegend(ax, position = :lt)

    configs = getParentState.(sectors,4)

    axConfs = [Axis(fig[2, i]; SW.getConfigAxis(configs[i])..., title = l,aspect = 1, xticklabelsvisible = false, yticklabelsvisible = false,xticks = 1:4,yticks = 1:4,titlecolor = colors[i]  ) for (i,l) in enumerate(labels)]

    for (ax,S) in zip(axConfs,configs)
        SW.plotSpinConfig!(ax,S)
    end
    rowsize!(fig.layout, 1, Relative(0.7))
    ylims!(ax, -0.34, 0)
    fig
end
##
function makeConf(UC,L,Spin)
    S = SW.stencilConfig(zeros(L,L),Spin;
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    S .= SW.getPeriodicState(UC,L,L)
    return S
end  
reducedConfigs = makeConf.(collect.(eachslice(SW.h5read("../../Data/reducedConfigs.h5","reducedConfigs"),dims=3)),4,1)

files =  [joinpath(root,file) for (root,_,files) in walkdir("/p/scratch/pmfrg/niggemann1/Spiderweb/DataS1_CT_RK_equil/SectorComp/L=16") for file in files]

mus = stack([h5read(file,"mu") for file in files])
perm = sortperm(mus)
mus = mus[perm]
ens = stack([h5read(file,"energies") for file in files])[:,perm]
Δens = stack([h5read(file,"Δenergies") for file in files])[:,perm]
##
with_theme(theme_SimpleTicks()) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"E_0 / N_\textrm{sites}")
    Nsites = 16^2
    for i in axes(ens,1)
        en = ens[i,:] ./ Nsites ./ (1 .- mu)
        deltaen = Δens[i,:]./ Nsites ./ (1 .- mu)

        scatterlines!(ax,mus,en;label = "Sector $i",marker = 'x')
        errorbars!(ax,mus,en,deltaen;label = "")
    end

    fig
end

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

with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 26,size = 1.5 .* (700,400))
    E_scal = 0.2
    Elin = E_scal .* (mus .-1)

    ax = Axis(fig[1,1],xlabel = L"\mu",ylabel = L"E_0/(N_\text{sites}(1 - \mu ))")
    Nsites = 12*12

 
    minsectors = argmin.(ens)

    text_locations = Point2f[]

    linewidths = LinRange(2,5.0, length(reducedConfigs))

    function getLS(width)
        if width < 4
            return :solid
        else
            return :dot
        end
    end
    # Generate all combinations of linestyles and linewidths

    for i_sector in eachindex(reducedConfigs)
        linewidth = linewidths[i_sector]
        linestyle = getLS(linewidth)
        
        en = ens[i_sector,:] ./ Nsites ./(1 .-mus)
        Δen = Δens[i_sector,:]./ Nsites ./(1 .-mus) 

        
        l = scatterlines!(ax,mus,en ,marker = '×';linestyle,linewidth,markersize = 0)
        errorbars!(ax,mus,en ,Δen,whiskerwidth=4)
        # text!(Point(-0.15,en[1]),text = "%$i_sector")

        # pos = round(Int,length(en)/length(AllSectors) * i_sector)
        # pos = max(1,pos)
        # pos = min(length(en),pos)

        # text!(ax,Point(mus[pos],en[pos]-0.01),text = L"%$i_sector", color = l.color[],fontsize = 10)
        # ppoint = isodd(i_sector) ? Point(mus[1]-0.05,en[begin]) : Point(mus[end]+0.05,en[end])
        ppoint = Point2f(minimum(mus)-0.04,en[begin])
        ppointRd = Point2f(ppoint[1],round(ppoint[2],digits = 2))
        
        while ppointRd in text_locations
            ppoint = Point2f(ppoint[1]+0.05,ppoint[2])
            ppointRd = Point2f(ppoint[1],round(ppoint[2],digits = 2))
        end
        push!(text_locations,ppointRd)
        
        text!(ax,ppoint,text = L"%$i_sector", color = l.color[],fontsize = 15,align = (:top,:top))
    end

    for i in eachindex(mus)[3:end]
        if minsectors[i-1] != minsectors[i]
            vlines!(ax,mus[i-1],color = :black,linestyle = :dash,linewidth = 0.8)
        end
    end
    colsize!(fig.layout,1,Relative(0.6))
    newfig = plotConfs!(fig[:,2:3],reducedConfigs,transpose=false)
    
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
