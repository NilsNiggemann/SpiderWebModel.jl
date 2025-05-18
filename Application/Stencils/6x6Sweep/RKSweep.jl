import SpiderWebModel as SW
using CairoMakie, MakieHelpers,Statistics
using SpiderWebModel.HDF5
using DataFrames
import SpiderWebModel.CircularArrays as CA
include("../plottingUtils.jl")
include("../FSSUtils.jl")

function photonDispersion(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
photonDispersion(k) = photonDispersion(k[1],k[2])

##
function pretty_scientific(x;kwargs...)
    if x == 0
        return L"0"
    end
    spl = split(string(x), "e")
    if length(spl) == 1
        return strd(x;kwargs...)
    end
    expo = floor(Int,log10(abs(x)))
    base = strd(x/(10. ^expo);kwargs...)
    return "$base×10^{$(expo)}"
end
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
res = getRes(ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/6x6Condensate_maxFlipInit_2/")
res2 = getRes(ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/6x6Condensate_maxFlipInit_3/")
filter!(x->!(x.mu ∈ res2.mu && x.L ∈ res2.L),res)
res = vcat(res,res2)
filter!(x->x.mu >=0.0,res)
sort!(res, [:mu, :L])

conf_6x6 = SW.periodicState6x6Condensate(30)
conf_6x6MF = copy(conf_6x6) .= h5read(ENV["MYSCRATCH"]*"Spiderweb/MaxFlip/6x6Condensate/L=30/Spin1GFMC_L=30_1.h5","maxConf")
# res = getRes(ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/S0*_longprop/")
##

with_theme(theme_SimpleTicks()) do 
    resNew = get_res(res2,mu=0.0,L=36)
    filter!(x->x.tau==0.2,resNew)
    println(resNew)
    ens = resNew.Energy
    tauax = only(unique(resNew.tau)) .* eachindex(ens[1])
    errlines(tauax[1:60],mean(ens)[1:60],std(ens)[1:60],axis = (;xlabel = L"τ",ylabel = L"E_0"))
    
end
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22, size = 200 .* (3, 2))
    ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"Energy")

    Linestyles = [:dash, :dot, :dashdot, :solid, :solid]
    scatterkwargs = Dict(
        16 => (;marker = '+'),
        20 => (;marker = '▴'),
        24 => (;marker = '●'),
        28 => (;marker = '×', markersize = 10),
        30 => (;marker = '▼', markersize = 10),
        36 => (;marker = '▲', markersize = 10),
        40 => (;marker = '■', markersize = 10),
    )

    for (L, linestyle) in zip((24,30,36), Linestyles)
        resL = get_res(res, L=L)
        mus = unique(resL.mu)
        energies = [getEnergy_tau(resL,8;mu) for mu in mus]

        enmean = mean.(energies) ./ L^2  #./ (1 .-mus)
        enstd = std.(energies) ./ L^2  #./ (1 .-mus)

        scatterlines!(ax, mus, enmean, label=L"L=%$L"; linestyle, scatterkwargs[L]...)
        errorbars!(ax, mus, enmean, enstd, whiskerwidth=5)
    end

    axislegend(ax, position=:lt)
    fig
end

##
with_theme(theme_SimpleTicks()) do
    L = 36
    resL = get_res(res,mu=0.0,L=L)
    filter!(x->x.tau==0.2,resL)
    SqsGFMC = resL.Sq

    SqMat = mean(SqsGFMC)
    SqErr = std(SqsGFMC)

    fig = Figure(size = 120 .* (4,4))
    ax = Axis(fig[1,1],xlabel = L"τ",ylabel = L"\mathcal{S}(\mathbf{q})")

    dTau = only(unique(resL.tau))
    tau = (axes(SqMat,3).-1) .*dTau
    # return heatmap(SqMat[:,:,20])
    Sq_examp = SqMat[:,:,10]
    inds = sort(collect(CartesianIndices(Sq_examp))[:],by = x->Sq_examp[x],rev=true)
    for I in ((2,1),(7,7),(10,3),(5,9),inds[1])
    # for I in inds[end-100:2:end-80]
    # for I in inds[[1,5,15,12,20,50,end÷4]]
    # for I in inds
        i,j = Tuple(I)
        range = 1:lastindex(tau)
        # scatterlines!(ax,tau[range],SqMat[i,j,range],marker = '×')
        # errorbars!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],whiskerwidth = 6,linewidth=0.5)
        errlines!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],linewidth=0.5)
    end
    fig
end

##
with_theme(theme_SimpleTicks()) do 
    L = 24
    muIndex = findfirst(>=(-0.2),res.mu)

    μ = res.mu[muIndex]

    SqsGFMC = SW.expand_Sq.(getSq(res,tau=12,mu=μ,L=L))
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
##
with_theme(theme_SimpleTicks()) do 
    L = 36
    muIndex = findfirst(>=(0.0),res.mu)

    μ = res.mu[muIndex]
    S_ref = SW.stencilConfig(zeros(L,L),1,boundaryCondition= :periodic) .= 2SW.periodicState6x6Condensate(L)
    SqFT_ref = SW.expand_Sq(abs2.(SW.FFTW.fft(S_ref)))

    SqsGFMC = SW.expand_Sq.(getSq(res,tau=14,mu=μ,L=L))
    SqMat = mean(SqsGFMC)
    SqErr = std(SqsGFMC)

    fig = Figure(size = 120 .* (4,4),fontsize = 22)

    xticks = yticks = PiTicks([0,pi])
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks)
    ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = false,yticklabelsvisible = false)

    Sq = SW.getSqCont(SqMat)

    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)

    heatmap!(ax,qx,qy,Sq)
    
    Sq_ref = SW.getSqCont(SqFT_ref)
    heatmap!(axFT,qx,qy,Sq_ref)

    fig
end

##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22, size = 200 .* (3, 2))
    ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"\Delta \mathcal{S}(\mathbf{q})")

    Linestyles = [:dash, :dot, :dashdot, :solid, :solid]
    scatterkwargs = Dict(
        16 => (;marker = '+'),
        20 => (;marker = '▴'),
        24 => (;marker = '●'),
        28 => (;marker = '×', markersize = 10),
        30 => (;marker = '▼', markersize = 10),
        36 => (;marker = '▲', markersize = 10),
        40 => (;marker = '■', markersize = 10),
    )

    for (L, linestyle) in zip((24,30,36), Linestyles)
        resL = get_res(res, L=L)
        mus = unique(resL.mu)
        norm_diffs = Float64[]
        norm_diff_errs = Float64[]
        qx = qy = trueMomenta(0, 2pi, L)
        for mu in mus
            SqsGFMC = SW.expand_Sq.(getSq(res, tau=14, mu=mu, L=L))
            SqMat = mean(SqsGFMC)
            SqErr = std(SqsGFMC)
            fittingCoefs = optimizeCoeffs(SqMat)
            SqFT = [SqFieldTheory(x, y, fittingCoefs...) for x in qx, y in qy]
            
            norm_diffs_individual = [SW.norm(Sq .- SqFT) ./length(Sq) for Sq in SqsGFMC]

            norm_diff = mean(norm_diffs_individual)
            norm_diff_err = std(norm_diffs_individual)
            push!(norm_diffs, norm_diff)
            push!(norm_diff_errs, norm_diff_err)
        end
        errorbars!(ax, mus, norm_diffs,norm_diff_errs, label=L"L=%$L",whiskerwidth = 10)
        scatterlines!(ax, mus, norm_diffs, label=L"L=%$L"; linestyle, scatterkwargs[L]...)
    end

    axislegend(ax, position=:rt,merge=true,unique=true)
    fig
end


##
resRandConfs = [getRes(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/RandConf_$i*_longprop") for i in 1:5]


function get_q_Cuts(L;numPoints=500)
    kpath = ["Γ","X","X'","Γ"]
    pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],numPoints)
    kpointlabels = Makie.latexstring.(kpath)
    tRange = eachindex(p1)
    
    qx = qy = trueMomenta(-0.5pi,1.5pi,L)

    xygrid = [(x,y) for x in qx, y in qy]

    xticks = (tRange[pointlabels],kpointlabels,)
    tRange_new,p1_discrete = rasterCurve(p1,xygrid,tRange)
    p1_points = xygrid[p1_discrete]
    return (;qx,qy,xygrid,p1_discrete,xticks,tRange_new,kpath,p1_points)
end
##
with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 22,size = 100 .*(3,3),  
    backgroundcolor = :transparent
    )
    ax_BCorr = Axis3(fig[1,1];xlabel = L"q_x",ylabel = L"q_y",zlabel = L"ω(\textbf{q})",
    zlabeloffset = 30,
    zticklabelpad = 10,
    xypanelcolor = :white,
    yzpanelcolor = :white,
    xzpanelcolor = :white,
    xypanelvisible = true,
    yzpanelvisible = true,
    xzpanelvisible = true,
    # xticks=PiTicks((-pi,0,pi)),yticks=PiTicks((-pi,0,pi)),
    xticks = PiTicks([0,pi]), yticks = PiTicks([0,pi]),
    zticks=SimpleTicks((0,1,2)),azimuth=1.6pi,elevation=0.3pi)
    qxPhot = trueMomenta(-0.5pi,1.5pi,250)
    qyPhot = trueMomenta(-0.5pi,1.5pi,250)
    # hm_B = surface!(ax_BCorr,qxPhot,qyPhot,photonDispersion,colormap = Makie.cgrad(:thermal,rev=false))
    hm_B = surface!(ax_BCorr,qxPhot,qyPhot,photonDispersion,colormap = Makie.cgrad(:inferno,rev=false))
    
    qx1,qy1 = qxPhot[1:15:end], qyPhot[1:15:end]
    # wireframe!(ax_BCorr,qx1,qy1,photonDispersion.(Iterators.product(qx1,qy1)),color = :black,linewidth = 0.5,overdraw = false)
    
    save("../../figs/PaperFigs/photonDispersion.png",fig)
    fig
    # text!(ax_BCorr,Point(pi/2,pi/2),text="TODO!",color = :black,align = (:center,:center),fontsize = 40)
end
##
function plotRandomCuts!(fig,resRandConfs)


    L_Plot = size(resRandConfs[1].Sq[1],1)

    (;qx,qy,xygrid,p1_discrete,xticks,tRange_new,kpath,p1_points) = get_q_Cuts(L_Plot;numPoints=100)

    axPath_Random = Axis(fig,ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = xticks,yticks = SimpleTicks()
    )

    combs = [
        (2,0.6),
        # (4,0.6),
        # (3,0.7),
        (4,0.8),
        (5,0.95),

    ]
    colors = Makie.Colors.distinguishable_colors(length(combs),parse.(Makie.Colors.RGB,[:black,:red,:blue,:orange]))
    

    labelpoints = [Point(0.5,0.15),Point(0.77,0.55),Point(0.85,0.8)]
    for (i,(Randsector,mu)) in enumerate(combs)
    # for (Randsector,mu) in Iterators.product(1:5,(0.7,0.8,0.9,0.95))
        SqsGFMC = SW.expand_Sq.(getSq(resRandConfs[Randsector],tau=20,mu=mu,L=36)) 
        SqMat = mean(SqsGFMC)
        SqErr = std(SqsGFMC)
        any(isnan,SqMat) && continue
        SqFunc = SW.getSqCont(SqMat)
        SqErrFunc = SW.getSqCont(SqErr)
        
        offset = 0.5*(i-1)
        Sqcut = [SqFunc(x,y) for (x,y) in p1_points] .+ offset
        Sqerrcut = [SqErrFunc(x,y) for (x,y) in p1_points]

        fittingCoefs = optimizeCoeffs(SqMat)
        
        A_fit,r_fit = pretty_scientific.(fittingCoefs,sigdigits=3)

        # text!(axPath_Random,labelpoints[i],text = L"A = %$(A_fit),\ r = %$(r_fit)",align = (:bottom,:right),space = :relative,color = colors[i],fontsize = 16)
        text!(axPath_Random,labelpoints[i],text = L"A = %$(A_fit)",align = (:bottom,:right),space = :relative,color = colors[i],fontsize = 16)
        text!(axPath_Random,labelpoints[i].-(0,0.1),text = L"r = %$(r_fit)",align = (:bottom,:right),space = :relative,color = colors[i],fontsize = 16)

        SqFT = [SqFieldTheory(q, fittingCoefs) for q in p1_points] .+ offset
        # return display(heatmap(qx,qy,(x,y) -> SqFieldTheory(SA[x,y],fittingCoefs)))
        l = lines!(axPath_Random, tRange_new, SqFT,label = L"\textrm{rand}_%$i\ μ = %$mu",color = colors[i])

        scatter!(axPath_Random, tRange_new, Sqcut,label =  L"\textrm{rand}_%$i\ μ = %$mu",color = colors[i])
        errorbars!(axPath_Random, tRange_new, Sqcut, Sqerrcut, whiskerwidth = 2, linewidth = 0.5,color = colors[i])
        
        # lines!(axPath_Random, [minimum(tRange_new),minimum(tRange_new)] .- 4offset, [0,offset], color = colors[i], linestyle = :dash)
        # lines!(axPath_Random, [maximum(tRange_new),maximum(tRange_new)].+ 4offset, [0,offset], color = colors[i], linestyle = :dash)
        hlines!(axPath_Random, [offset], [minimum(tRange_new),maximum(tRange_new)], color = (colors[i],0.3), linestyle = :solid,linewidth = 0.8,)
    end
    # axislegend(axPath_Random,position = :lt,merge=true,nbanks=3,labelsize=18)
    return axPath_Random
end



using CairoMakie.FileIO
with_theme(theme_SimpleTicks()) do 

    fig = Figure(fontsize = 22,size = 550 .*(3,1))
    
    fig[1:2,1] = ConfPanels  = GridLayout()
    fig[1:2,2:4] = SqPanels  = GridLayout()
    fig[1:2,5] = FTPanels = GridLayout()
    fig[1,6:8] = ThirdCol = GridLayout()
    fig[2,6:8] = FourthCol = GridLayout()
    
    # boxSq = Box(fig[1:2, 1:3], linestyle = :solid,alignmode = Mixed(left = -70, right = -12, top = -8, bottom = -35),color = (:grey,0.1),strokecolor = :black,cornerradius = 10)
    # boxFT = Box(fig[1:2, 4], linestyle = :solid,alignmode = Mixed(left = -40, right = -12, top = -8, bottom = -35),color = (:grey,0.1),strokecolor = :black,cornerradius = 10)

    # boxScal = Box(fig[3, 1:4], linestyle = :solid,alignmode = Mixed(left = -70, right = -12, top = -12, bottom = -60),color = (:grey,0.1),strokecolor = :black,cornerradius = 10)
    # boxRand = Box(fig[4, 1:4], linestyle = :solid,alignmode = Mixed(left = -70, right = -12, top = -35, bottom = -60),color = (:grey,0.1),strokecolor = :black,cornerradius = 10)

    # Makie.translate!(boxSq.blockscene, 0, 0, -100)
    # Makie.translate!(boxFT.blockscene, 0, 0, -100)
    # Makie.translate!(boxScal.blockscene, 0, 0, -100)
    # Makie.translate!(boxRand.blockscene, 0, 0, -100)
    # Makie.translate!(box2.blockscene, 0, 0, -100)

    SqPanels[1,1:3] = Sq_Heatmaps = GridLayout()
    FTPanels[1,1] = FT_DispFig = GridLayout()
    SqPanels[2,1:3] = SqCuts = GridLayout()
    FTPanels[2,1] = FT_SqFig = GridLayout()
    ThirdCol[1,1] = FSS_Plot = GridLayout()
    FourthCol[1,1] = SqRandomCuts = GridLayout()


    ax_conf1 = Axis(ConfPanels[1,1];SW.getConfigAxis(conf_6x6)...,xticks = [10,20,30],yticks = [10,20,30],xminorgridwidth = 0.1, xticklabelsvisible=false,yminorgridwidth = 0.1)
    ax_conf2 = Axis(ConfPanels[2,1];SW.getConfigAxis(conf_6x6MF)...,xticks = [10,20,30],yticks = [10,20,30],xminorgridwidth = 0.1, yminorgridwidth = 0.1)
    linkxaxes!(ax_conf1,ax_conf2)
    SW.plotSpinConfig!(ax_conf1,conf_6x6,constraintkwargs = (;markersize = 7))
    SW.plotSpinConfig!(ax_conf2,conf_6x6MF,constraintkwargs = (;markersize = 7))
    PiTicksArgs = (;xticks = PiTicks([0,pi]), yticks = PiTicks([0,pi]))

    L_Plot = 36

    (;qx,qy,xygrid,p1_discrete,xticks,tRange_new,kpath,p1_points) = get_q_Cuts(L_Plot)

    FSS_Plot[1,1] = ax_scal = Axis(fig,xlabel = L"μ",ylabel = L"\textrm{max}(\mathcal{S}(\mathbf{q}))")
    # FSS_Plot[1,1] = ax_scal = Axis(fig,xlabel = L"μ",ylabel = L"\textrm{max}_\Delta \mathcal{S}(\mathbf{q})")
    FSS_Plot[1,1] = ax_scal2 = Axis(fig,yaxisposition = :right,ylabel = L"$\delta \mathcal{S}(\mathbf{q})$ (dashed)",yticklabelcolor = :gray20,ylabelcolor = :gray20,xgridvisible = false,ygridvisible = false)
    # ylims!(ax_scal,1.15,1.4)
    # linkxaxes!(ax_scal,ax_scal2)
    # mu_show = (0.8,0.8,0.8)
    mu_show = (0.0,0.4,0.8)
    # mu_show = (-0.2,0.0,0.2)

    SqsGFMC = [SW.expand_Sq.(getSq(res,tau=7,mu=mu,L=L_Plot)) for mu in mu_show]
    SqMat = mean.(SqsGFMC)
    SqErr = std.(SqsGFMC)
    
    framecolors = [:black,:grey,:magenta]
    # framecolors = [:black,:blue,:red]
    spinewidth= 3
    with_theme(theme_PiTicks()) do 
        Sq_Heatmaps[1,1] = ax_mu1 = Axis(fig;xlabel = L"q_x",ylabel = L"q_y",aspect=1,PiTicksArgs...,xlabelpadding = -10,spinewidth,spinecolors(framecolors[1])...,)
        Sq_Heatmaps[1,2] = ax_mu2 = Axis(fig;xlabel = L"q_x",aspect=1,yticklabelsvisible = false,PiTicksArgs...,xlabelpadding = -10,spinewidth,spinecolors(framecolors[2])...,)
        Sq_Heatmaps[1,3] = ax_mu3 = Axis(fig;xlabel = L"q_x",aspect=1,yticklabelsvisible = false,PiTicksArgs...,xlabelpadding = -10,spinewidth,spinecolors(framecolors[3])...,)
        
        # ax_mu3.alignmode = Mixed(bottom = 30)
        linkyaxes!(ax_mu1,ax_mu2,ax_mu3)


        colorrange = extrema(stack(SqMat))

        for (i,ax) in enumerate((ax_mu1,ax_mu2,ax_mu3))
            SqCont = SW.getSqCont(SqMat[i])
            hm = heatmap!(ax,qx,qy,SqCont,colormap = :viridis;
            colorrange
            )
        end

        fittingCoefs = optimizeCoeffs(SqMat[end])
        SqFT(qx,qy) = SqFieldTheory(qx, qy, fittingCoefs...)

        # A_fit,r_fit = pretty_scientific.(fittingCoefs,sigdigits = 3)
        # rinv = strd(1/fittingCoefs[2],sigdigits = 3)
        # rinv = pretty_scientific(1/fittingCoefs[2])
        # ax_FT = Axis(FSS_Plot[1,1],width=Relative(0.3),height =Relative(0.3),halign = 0.9,valign = 0.9;aspect=1,title = L"$A = %$(A_fit),\ r = %$(r_fit)$",xlabel = L"q_x",ylabel = L"q_y",PiTicksArgs...,xlabelpadding=-5.)
        # ax_FT = Axis(FT_SqFig[1,1];aspect=1,title = L"$A = %$(A_fit),\ r = %$(r_fit)$",xlabel = L"q_x",ylabel = L"q_y",PiTicksArgs...,xlabelpadding=-5.)
        # ax_FT = Axis(FT_SqFig[1,1];aspect=1,title = L"$A = %$(A_fit),\ r^{-1} = %$(rinv)$",xlabel = L"q_x",ylabel = L"q_y",PiTicksArgs...,xlabelpadding=-5.,spinewidth,spinecolors(framecolors[3])...)
        ax_FT = Axis(FT_SqFig[1,1];aspect=1,xlabel = L"q_x",ylabel = L"q_y",PiTicksArgs...,xlabelpadding=-8,ylabelpadding=-3,spinewidth,spinecolors(framecolors[3])...)
        hmFT = heatmap!(ax_FT, qx, qy, SqFT;colorrange)

        linepoints = getindex.(Ref(KPoints), kpath)
        lines!(ax_FT, linepoints, color = :red,linestyle = :dash, linewidth = 1.5)
        for label in kpath
            Q = KPoints[label]
            text!(ax_FT, Point(Q...), text=Makie.latexstring(label), strokecolor=(:black,0.3), align=(:center, :center),strokewidth=4)
            text!(ax_FT, Point(Q...), text=Makie.latexstring(label), color=:white, align=(:center, :center))
        end
        for (i, ax) in enumerate((ax_mu1, ax_mu2, ax_mu3,))
            text!(ax, Point(0.5,0.99), text=L"\mu = %$(mu_show[i])", align=(:center, :top),space = :relative, strokecolor = (:black,0.3), strokewidth = 4)
            text!(ax, Point(0.5,0.99), text=L"\mu = %$(mu_show[i])", align=(:center, :top),space = :relative,color = :white)
        end
        cb = Colorbar(Sq_Heatmaps[0,1:3],ticks = SimpleTicks(),width = Relative(1),height = 10;colorrange,vertical=false,label = L"\mathcal{S}(\mathbf{q})",labelpadding = -2)
        cb.alignmode = Mixed(top = 15)
        # cb.alignmode = Mixed(right = 0,top = -40,bottom = -20)
    end
    with_theme(theme_SimpleTicks()) do 
        SqCuts[1,1] = axPath1 = Axis(fig;
        # ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" ,
         xticks = xticks,
         spinewidth,spinecolors(framecolors[1])...,
         backgroundcolor = :white,
         ylabel = L"\mathcal{S}(\mathbf{q})",
        )
        SqCuts[1,2] = axPath2 = Axis(fig;
        # xlabel = L"\mathbf{q}",
        spinewidth,spinecolors(framecolors[2])...,
         xticks = xticks,yticklabelsvisible = false,
         backgroundcolor = :white,
        )
        SqCuts[1,3] = axPath3 = Axis(fig;
        # xlabel = L"\mathbf{q}",
        spinewidth,spinecolors(framecolors[3])...,
         xticks = xticks,yticklabelsvisible = false,
         backgroundcolor = :white,
        )
        linkyaxes!(axPath1,axPath2,axPath3)

        for (i,ax,SqMat,SqErr) in zip(eachindex(SqsGFMC), (axPath1, axPath2, axPath3), SqMat, SqErr)
            SqFunc = SW.getSqCont(SqMat)
            SqErrFunc = SW.getSqCont(SqErr)
            Sqcut = [SqFunc(x,y) for (x,y) in p1_points]
            Sqerrcut = [SqErrFunc(x,y) for (x,y) in p1_points]

            fittingCoefs = optimizeCoeffs(SqMat)

            begin
                A,r = pretty_scientific.(fittingCoefs,sigdigits = 3)
                # rinv = strd(1/fittingCoefs[2],sigdigits = 3)
                # rinv = pretty_scientific(1/fittingCoefs[2])
                text!(ax, Point(0.95,0.98), text = L"A = %$(A)", align = (:right,:top),space = :relative,color = framecolors[i],fontsize = 16)
                # text!(ax, Point(0.95,0.88), text = L"r^{-1} = %$(rinv)", align = (:right,:top),space = :relative,color = framecolors[i],fontsize = 16)
                text!(ax, Point(0.95,0.88), text = L"r = %$(r)", align = (:right,:top),space = :relative,color = framecolors[i],fontsize = 16)
            end


            SqFT = [SqFieldTheory(q, fittingCoefs) for q in p1_points]
            lines!(ax, tRange_new, SqFT, color = :red, linestyle = :dash)
            scatter!(ax, tRange_new, Sqcut, marker = :circle, markersize = 5, color = :black)
            errorbars!(ax, tRange_new, Sqcut, Sqerrcut, whiskerwidth = 2, linewidth = 0.5, color = :black)
            ylims!(ax,0,maximum(Sqcut) + 0.5)

        end
    end

    axPathsRandom = SqRandomCuts[1,1] = plotRandomCuts!(fig, resRandConfs)
    Legend(SqRandomCuts[1,:,Top()],axPathsRandom,merge=true,nbanks=3,padding= (10,10,2,2),tellheight=false,margin = (0,0,60,30),labelsize=22)

    with_theme(theme_PiTicks()) do 
        FT_DispFig[1,1] = ax_disp = Axis(fig,
        ylabelvisible = false,
        aspect = DataAspect(),
        xlabelvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksvisible = false,
        backgroundcolor = (:white, 0.0),
        # backgroundcolor = (:black, 0.8),
        yticksvisible = false,
        xgridvisible = false,
        tellheight = false,
        tellwidth = false,    
        ygridvisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        )
        img = load("../../figs/PaperFigs/photonDispersion.png")[90:end,105:end]
        image!(ax_disp, rotr90(img),transparency=true)
        ax_disp.alignmode = Mixed(left = -35,right = -15, top = -0,bottom = -10)
    end
    Linestyles = [:solid, :solid, :solid]
    scatterkwargs = Dict(
        16 => (;marker = '+'),
        20 => (;marker = '▴'),
        24 => (;marker = '●'),
        28 => (;marker = '×', markersize = 10),
        30 => (;marker = '▼', markersize = 10),
        36 => (;marker = '▲', markersize = 10),
        40 => (;marker = '■', markersize = 10),
    )
    ylim_max = 0

    cols = Dict(
        24 => :blue,
        30 => :red,
        36 => :black,
    )
    for (L, linestyle) in zip((24,30,36), Linestyles)
        resL = get_res(res, L=L)
        isempty(resL) && continue
        mus = unique(resL.mu)
        norm_diffs = Float64[]
        norm_diff_errs = Float64[]
        max_Sqs = Float64[]
        max_Sqs_err = Float64[]
        
        qx = qy = trueMomenta(0, 2pi, L)
        for mu in mus
            taufilter = mu ==0 ? 1 : 11
            SqsGFMC = SW.expand_Sq.(getSq(res, tau=taufilter, mu=mu, L=L))
            SqMat = mean(SqsGFMC)
            SqErr = std(SqsGFMC)
            fittingCoefs = optimizeCoeffs(SqMat)

            SqFT = [SqFieldTheory(x, y, fittingCoefs...) for x in qx, y in qy]
            
            norm_diffs_individual = [SW.norm(Sq .- SqFT) ./L^2 for Sq in SqsGFMC]

            # max_deviation_individual = [maximum(abs.(Sq .- SqFT)) for Sq in SqsGFMC]
            max_Sq_individual = [maximum(Sq) for Sq in SqsGFMC]
            # max_Sq_individual = [getxi(CA.CircularArray(Sq)) for Sq in SqsGFMC]
            
            norm_diff = mean(norm_diffs_individual)
            norm_diff_err = std(norm_diffs_individual)
            
            push!(norm_diffs, norm_diff)
            push!(norm_diff_errs, norm_diff_err)

            max_Sq = mean(max_Sq_individual)
            max_Sq_err = std(max_Sq_individual)
            push!(max_Sqs, max_Sq)
            push!(max_Sqs_err, max_Sq_err)
        end

        begin
            errorbars!(ax_scal2, mus, norm_diffs,norm_diff_errs, label=L"L=%$L",whiskerwidth = 10,color = cols[L])
            scatterlines!(ax_scal2, mus, norm_diffs, label=L"L=%$L"; linestyle=:dash,color = cols[L], linewidth = 0.5,scatterkwargs[L]...)
            # scatter!(ax_scal2, mus, norm_diffs, label=L"L=%$L";color = cols[L], scatterkwargs[L]...)
        end
        # ylim_max = max(ylim_max,maximum(max_Sqs))
        errorbars!(ax_scal, mus, max_Sqs,max_Sqs_err, label=L"L=%$L",whiskerwidth = 10,color = cols[L])
        scatterlines!(ax_scal, mus, max_Sqs, label=L"L=%$L"; linestyle,color = cols[L], scatterkwargs[L]...)
    end
    
    # ylims!(ax_scal,nothing,1.6ylim_max)
    axislegend(ax_scal,position = :rt,merge=true)
    
    rowsize!(SqPanels,1,Relative(0.6))
    rowsize!(SqPanels,2,Relative(0.4))
    rowsize!(fig.layout,1,Relative(0.6))
    rowsize!(fig.layout,2,Relative(0.4))
    colsize!(fig.layout,1,Relative(0.2))
    # rowsize!(fig.layout,1,Relative(0.22))
    # rowsize!(fig.layout,2,Relative(0.37))
    # rowsize!(fig.layout,3,Relative(0.22))
    # rowsize!(fig.layout,4,Relative(0.2))
    # rowsize!(SqCuts,1,Relative(0.8))
    # rowgap!(Sq_Heatmaps,1,5)
    # rowgap!(SqPanels,1,-10)

    # colsize!(fig.layout,4,Relative(0.4))

    # rowsize!(SqPanels,1,Relative(0.6))
    # rowsize!(SqPanels,2,Relative(0.5))
    # rowsize!(FTPanels,1,Relative(0.6))
    # rowsize!(FTPanels,2,Relative(0.34))
    # colgap!(Sq_Heatmaps,1,5)
    # colgap!(Sq_Heatmaps,2,5)
    # colgap!(SqCuts,1,5)
    # colgap!(SqCuts,2,5)
    # rowgap!(fig.layout,1,-300)
    # rowgap!(fig.layout,2,10)
    # rowgap!(fig.layout,3,50)

    ax_scal.alignmode = Mixed(right = 87)

    ax_scal2.alignmode = Mixed(right = 0)
    # Label(fig[1,1,TopLeft()],L"(a)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,30,-20,0))
    # Label(SqCuts[1,1,TopLeft()],L"(b)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,30,0,0))
    # Label(fig[1,4,TopLeft()],L"(c)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,0,-20,0))
    # Label(FT_SqFig[1,1,TopLeft()],L"(d)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,0,0,0))
    # Label(fig[3,1,TopLeft()],L"(e)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,30,-10,0))
    # Label(fig[4,1,TopLeft()],L"(f)$$", fontsize = 24,tellheight=false,tellwidth=false,padding = (0,30,30,0))
    save("../../figs/SqFieldTheoryComparison.pdf",fig)
    fig
    
end
