import SpiderWebModel as SW
using CairoMakie, MakieHelpers,Statistics
using SpiderWebModel.HDF5
using DataFrames

include("../plottingUtils.jl")

##
function file_format(filename)
    allkeys_1 = ["energies","mu","tau","SqsGFMC","p_Sq"]
    allkeys_2 = ["Energy","mu","τ","StructureFactor"]
    h5open(filename,"r") do file
        all(k->haskey(file,k),allkeys_1) && return 1
        all(k->haskey(file,k),allkeys_2) && return 2
        return 0
    end
end

function getSq_tau(res::DataFrame,tau)
    Dtaus = res.tau

    taus = [axes(Sqs,3) .* t for (Sqs,t) in zip(res.Sq,Dtaus)]
    tauInds = [findfirst(>=(tau),t) for t in taus]

    # return tauInds
    return [res.Sq[i][:,:,ti] for (i,ti) in enumerate(tauInds)]
end

function getEnergy_tau(res::DataFrame,tau;kwargs...)
    res = get_res(res;kwargs...)
    Dtaus = res.tau

    taus = [axes(Energies,1) .* t for (Energies,t) in zip(res.Energy,Dtaus)]
    tauInds = [findfirst(>=(tau),t) for t in taus]

    # return tauInds
    return [res.Energy[i][ti] for (i,ti) in enumerate(tauInds)]
end

function get_res(res::DataFrame;mu=nothing,L=nothing)
    mu_func = isnothing(mu) ? x->true : x->x.mu == mu
    L_func = isnothing(L) ? x->true : x->x.L == L
    res_new = filter(row -> mu_func(row) && L_func(row),res)
end

function getSq_tau(res::DataFrame,::Nothing)
    return res.Sq
end

function getSq(res::DataFrame;tau=nothing,mu=nothing,L=nothing)
    res_new = get_res(res,mu=mu,L=L)
    return getSq_tau(res_new,tau)
end

function getSq_tau_mean_std(res,tau)
    Sqtau = getSq_tau(res,tau)
    dropmean(Sqtau,dims=4), dropstd(Sqtau,dims=4)
end

function getSq_tau_mean_std(res::DataFrame,tau)
    Sqtau = getSq_tau(res,tau)
    mean(Sqtau), std(Sqtau)
end

function getRes(folder)
    files = let
        filesunsrt = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
        fileformats = file_format.(filesunsrt)
        invalid_files = findall(iszero,fileformats)
        if !isempty(invalid_files)
            println("invalid files:")
            println(filesunsrt[invalid_files])
        end
        filesunsrt = [f for (f,i) in zip(filesunsrt,fileformats) if i == 2]
        mus = [h5read(file,"mu") for file in filesunsrt]
        filesunsrt[sortperm(mus)]
    end

    L = [h5read(file,"L") for file in files]
    mu = [h5read(file,"mu") for file in files]
    Energy = [h5read(file,"Energy") for file in files]
    Sq = [h5read(file,"StructureFactor") for file in files]
    tau = [h5read(file,"τ") for file in files]

    res= DataFrame(;L,mu,Energy,Sq,tau,files)
    sort!(res,[:mu,:L])
end
##
resRandConfs = [getRes(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_$i*_longprop") for i in 1:5]
##
with_theme(theme_SimpleTicks()) do 
    resNew = get_res(resRandConfs[2],mu=0.6,L=36)

    ens = resNew.Energy
    tauax = only(unique(resNew.tau)) .* eachindex(ens[1])
    errlines(tauax,mean(ens),std(ens),axis = (;xlabel = L"τ",ylabel = L"E_0"))
    
end
##

function plotSqOverview(res,mu)
    with_theme(theme_SimpleTicks()) do 
        L = 36
        muIndex = findfirst(>=(mu),res.mu)

        μ = res.mu[muIndex]

        SqsGFMC = SW.expand_Sq.(getSq(res,tau=20,mu=μ,L=L))
        SqMat = mean(SqsGFMC)
        SqErr = std(SqsGFMC)
        fittingCoefs = optimizeCoeffs(SqMat)

        fig = Figure(size = 120 .* (4,4),fontsize = 22)

        xticks = yticks = PiTicks([0,pi])
        axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks)

        ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = false,yticklabelsvisible = false)

        # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
        Sq = SW.getSqCont(SqMat)
        Sqerr = SW.getSqCont(SqErr)
        qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
        heatmap!(ax,qx,qy,Sq)

        # SqFT = [SqFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]
        heatmap!(axFT,qx,qy,(qx,qy)->SqFieldTheory(qx,qy,fittingCoefs...))

        colors = (:red,:blue,:magenta)


        colorFT = :black
        colorGFMC = :red

        kpath = ["Γ","X","X'","Γ"]
        pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],100)
        kpointlabels = Makie.latexstring.(kpath)
        tRange = eachindex(p1)
        xygrid = [(x,y) for x in qx, y in qy]

        axPath = Axis(fig[2,1:2],ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
        )
        tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)


        p1_points = xygrid[p1_discrete]

        Sqcut = [Sq(x,y) for (x,y) in p1_points]
        Sqerrcut = [Sqerr(x,y) for (x,y) in p1_points]
        SqFT = [SqFieldTheory(q,fittingCoefs) for q in p1_points]

        # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
        scatter!(ax,p1_points,marker = '∘' ,color = colorGFMC,markersize = 15)
        scatterlines!(axFT,p1_points,color = colorFT,linestyle = :dash,marker = '●',markersize = 2)
        # tRange = SW.norm.(p1).^2
        scatterlines!(axPath,tRange,SqFT,color = colorFT,linestyle = :dash,marker = '●',markersize = 8)

        text!(axFT,Point(0,0),text="Γ",color = :white,align = (:center,:center))
        text!(axFT,Point(pi,0),text="X",color = :white,align = (:center,:center))
        text!(axFT,Point(0,pi),text="X'",color = :white,align = (:center,:center))

        scatter!(axPath,tRange,Sqcut,
        marker = '∘',markersize = 18,color = colorGFMC)
        errorbars!(axPath,tRange,Sqcut,Sqerrcut,color = colorGFMC,whiskerwidth = 6,linewidth=0.5)

        rowsize!(fig.layout,1,Relative(0.5))
        text!(axPath,Point(tRange[end-end÷8],1.2),text=L"μ = %$μ",align = (:center,:center))
        fig
    end
end

##
display.(plotSqOverview.(resRandConfs,0.6));
##
display.(plotSqOverview.(resRandConfs,0.8));
##
display.(plotSqOverview.(resRandConfs,0.95));
##
plotSqOverview(resRandConfs[5],0.6)

##
plotSqOverview(resRandConfs[4],0.8)