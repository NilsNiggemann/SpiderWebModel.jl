import SpiderWebModel as SW
using CairoMakie, MakieHelpers,Statistics
using SpiderWebModel.HDF5
using DataFrames

include("../plottingUtils.jl")
include("../FSSUtils.jl")

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

function getRes_2(folder)
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
res = getRes_2("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/StairCase/L=28/")
##
with_theme(theme_SimpleTicks()) do
    ind = 1
    L = 28
    resL = get_res(res,L=L,mu=0.8)

    Nsites = length(getSq_tau(res,10.))
    enmean = mean(resL.Energy)./ Nsites

    enstd = std(resL.Energy)./ Nsites
    # enmean = mean(entest,dims=2)[:,1]
    # enstd = std(entest,dims=2)[:,1]
    tau = only(unique(resL.tau)) .* (eachindex(enmean) .-1)
    errlines(tau,enmean,enstd,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # lines((eachindex(enmean).-1) .*tau,enmean,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # band!((eachindex(enmean).-1) .*tau,enmean - enstd , enmean + enstd,color = (:black,0.2))
    # current_figure()

end
##
#
with_theme(theme_PiTicks()) do
    L = 28
    muPlot = [0.1,0.8,0.9,1.0]
    
    fig = Figure(fontsize = 22,size = 170 .*(length(muPlot),1.85))
    ticks = PiTicks([0,pi])
    axes = [Axis(fig[1,i],aspect=1,
    # title = L"μ = %$(muPlot[i])",
    yticklabelsvisible=i==1,xticks=ticks,yticks=ticks,xlabel = L"q_x",ylabel = L"q_y", ylabelvisible = i==1) for i in eachindex(muPlot)]
    
    for (i,mu) in enumerate(muPlot)
        # return getSq(res,tau=15.,L=L,mu=mu)
        SqMat = SW.expand_Sq(mean(getSq(res,tau=15.,L=L,mu=mu)))

        kx = ky = trueMomenta(-pi/2,1.5pi,L)

        SqFunc = SW.getSqCont(SqMat)
        Sqpl = SqFunc.(Iterators.product(kx,ky))
        hm = heatmap!(axes[i],kx,ky,Sqpl,colormap = :viridis)
        Colorbar(fig[0,i],hm,ticks = SimpleTicks(),vertical = false,flipaxis = true,label = L"\mathcal{S}(\mathbf{q})",width = Relative(0.9))
        
        band!(axes[i],[-0.5pi,0.5pi],[1.1pi,1.1pi],[1.5pi,1.5pi],color = (:black,0.5))
        text!(axes[i],Point(0,1.3pi),text=L"μ = %$mu",color = :white,align = (:center,:center))
        
        # text!(axes[i],Point(0,1.3pi),text=L"μ = %$mu",color = :black,align = (:center,:center))
    end
    rowgap!(fig.layout,1,5)
    colgap!(fig.layout,1,0)
    colgap!(fig.layout,2,0)
    colgap!(fig.layout,3,0)

    fig
end
##
with_theme(theme_SimpleTicks()) do
    L = 28
    muIndex = findfirst(>=(0.2),res[L].mus)
    SqsGFMC = res[L].Sqs[:,:,:,:,muIndex]./ 4
    SqMat = dropmean(SqsGFMC,dims=4)
    SqErr = dropstd(SqsGFMC,dims=4)
    fig = Figure(size = 120 .* (4,4))
    ax = Axis(fig[1,1],xlabel = L"τ",ylabel = L"\mathcal{S}(\mathbf{q})")
    p_Sq = res[L].p_Sq[:,muIndex]
    dTau = res[L].taus[muIndex]
    tau = p_Sq .*dTau
    # return heatmap(SqMat[:,:,20])
    Sq_examp = SqMat[:,:,10]
    inds = sort(collect(CartesianIndices(Sq_examp))[:],by = x->Sq_examp[x],rev=true)
    # for I in ((5,5),(7,7),(10,3),(5,9))
    for I in inds[[1,5,15,12,20,50]]
        i,j = Tuple(I)
        range = 1:120
        # scatterlines!(ax,tau[range],SqMat[i,j,range],marker = '×')
        # errorbars!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],whiskerwidth = 6,linewidth=0.5)
        errlines!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],linewidth=0.5)
    end
    fig
end
##
with_theme(theme_SimpleTicks()) do 
    L = 24
    muIndex = findfirst(>=(0.7),res[L].mus)
    SqsGFMC = res[L].Sqs[:,:,100,:,muIndex]./ 4
    SqMat = dropmean(SqsGFMC,dims=3)
    SqErr = dropstd(SqsGFMC,dims=3)
    fittingCoefs = optimizeCoeffs(SqMat)

    μ = res[L].mus[muIndex]
    fig = Figure(size = 120 .* (4,4))

    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks(), yticks = PiTicks())

    ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks(), yticks = PiTicks(),ylabelvisible = false,yticklabelsvisible = false)

    # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
    Sq = SW.getSqCont(SqMat)
    Sqerr = SW.getSqCont(SqErr)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    
    SqFT = [SqFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]
    heatmap!(axFT,qx,qy,SqFT)

    q_path(r,phi) = (r*cos(phi),r*sin(phi))
    qr = LinRange(0,.35pi,100)
    
    colors = (:red,:blue,:magenta)
    


    let color = :black

        KPoints = Dict([
            "Γ" => SVector(0,0),
            "X" => SVector(pi,0),
            "M" => SVector(pi,pi),
            "X'" => SVector(0,pi)
            ])
        kpath = ["Γ","X","X'","Γ"]
        pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],500)
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
        scatter!(ax,p1_points,marker = '∘' ,color = color,markersize = 10)
        scatterlines!(axFT,p1_points,color = color,linestyle = :dash,marker = '●',markersize = 2)
        # tRange = SW.norm.(p1).^2
        scatter!(axPath,tRange,Sqcut,
        marker = '∘',markersize = 15,color = color)
        errorbars!(axPath,tRange,Sqcut,Sqerrcut,color = color,whiskerwidth = 6,linewidth=0.5)
        scatterlines!(axPath,tRange,SqFT,color = color,linestyle = :dash,marker = '●',markersize = 4)
        
        text!(axFT,Point(0,0),text="Γ",color = :white,align = (:center,:center))
        text!(axFT,Point(pi,0),text="X",color = :white,align = (:center,:center))
        text!(axFT,Point(0,pi),text="X'",color = :white,align = (:center,:center))
    end

    # for (phi,color) in zip([0,pi/4],colors)
    #     qpoints_raw = q_path.(qr,phi)
    #     qpoints = sort!(unique!(roundToTrueMomenta.(qpoints_raw,size(SqMat,1)-1)), by = SW.norm)

    #     Sqcut = Sq.(qpoints)
    #     Sqerrcut = Sqerr.(qpoints)
        
    #     # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
    #     SqFT = [SqFieldTheory(q,fittingCoefs...) for q in qpoints]
    #     scatter!(ax,qpoints,marker = '×' ,color = color)
    #     scatterlines!(axFT,Point.(qpoints),color = color,linestyle = :dash,marker = '●',markersize = 4)
    #     qnorms_sq = SW.norm.(qpoints).^2
    #     scatter!(ax2,qnorms_sq,Sqcut,
    #     marker = '×',markersize = 15,color = color)
    #     errorbars!(ax2,qnorms_sq,Sqcut,Sqerrcut,color = color,whiskerwidth = 6,linewidth=0.5)
    #     scatterlines!(ax2,qnorms_sq,SqFT,color = color,linestyle = :dash,marker = '●',markersize = 4)
    # end
    rowsize!(fig.layout,1,Relative(0.5))

    Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end

with_theme(theme_PiTicks()) do 
    mu = 0.7
    fig = Figure(fontsize = 22,size = 400 .*(1.4,2.4))
    ticks = PiTicks([0,pi])
    axes = [
        Axis(fig[1,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = %$(mu)"),
        Axis(fig[2,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = %$(mu)"),
    ]

    for (i,ax1,L) in zip((1,2),axes,(20,24))
        # Sq = getSq_tau(res[L],13)
        mus = res[L].mus
        
        muIndex = findfirst(>=(mu),mus)
        SqsGFMC = res[L].Sqs[:,:,end,:,muIndex]./ 4
        # SqsGFMC = Sq[:,:,:,muIndex]./ 4

        SqMat = dropdims(mean(SqsGFMC,dims=3),dims=3)
        SqErr = dropdims(std(SqsGFMC,dims=3),dims=3)


        # ax2 = Axis(fig[1,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = 0.6")
        kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
        SqFunc = SW.getSqCont(SqMat)

        # SqFunc = SW.getSqCont(SqMat)

        # hm = heatmap!(ax1,kx,ky,SqFunc.(Iterators.product(kx,ky)),colormap = :viridis)
        hm = heatmap!(ax1,kx,ky,SqFunc.(Iterators.product(kx,ky)),colormap = :viridis)
        xi = getxi(SqMat)
        @info "" L maximum(SqMat) /L^2 xi xi/L mus[muIndex]
        Colorbar(fig[i,2],hm)
    end
    fig
end

##
function getXiLs(res,res_36)
    
    xis = Dict{Int,Vector{Float64}}()
    xis_err = Dict{Int,Vector{Float64}}()
    mus = Dict{Int,Vector{Float64}}()

    let L=36
        unique_mus = unique(res_36.mu)
        k = trueMomenta(0,2pi,L)
        i_k = findfirst(==(pi/2),k)

        Sq = [getSq(res_36,mu = mu,L = L,tau=10) for mu in unique_mus]
        xi_L = [getXis(Sq,CartesianIndex(i_k,i_k)) for Sq in Sq]
        xis[L] = mean.(xi_L) ./ L
        xis_err[L] = std.(xi_L) ./ L
        mus[L] = unique_mus
    end

    for L in (20,24,28)
        Sqs = eachslice(getSq_tau(res[L],10),dims=(3,4))
        k = trueMomenta(0,2pi,L)
        i_k = findfirst(==(pi/2),k)
        xi_L = [getXis(Sq,CartesianIndex(i_k,i_k)) for Sq in eachslice(Sqs,dims=1)]
        xis[L] = mean(xi_L) ./ L
        xis_err[L] = std(xi_L) ./ L
        mus[L] = res[L].mus
    end

    return (;xis,xis_err,mus)
end
xi_res = getXiLs(res,res_36)



xis_intPol = let

    Dict(L=> interpolate((xi_res.mus[L],), xi_res.xis[L], Gridded(Linear())) for L in keys(xi_res.xis))
end
crossings = detect_crossings(xis_intPol)
##
with_theme(theme_SimpleTicks()) do 

    fig = Figure(fontsize = 22,size = 400 .*(3,2))
    
    FSS_Plot = GridLayout()
    Sq_Heatmaps = GridLayout()
    SqCuts = GridLayout()
    SqRandomCuts = GridLayout()
    BCorrPlot = GridLayout()

    fig.layout[1:2,1:3] = FSS_Plot
    fig.layout[1,4:7] = Sq_Heatmaps
    fig.layout[2,4:7] = SqCuts
    fig.layout[3,1:7-1] = SqRandomCuts
    fig.layout[3,7] = BCorrPlot
    
    PiTicksArgs = (;xticks = PiTicks([0,pi]), yticks = PiTicks([0,pi]))

    
    L_Plot = 36
    qx = qy = trueMomenta(-0.5pi,1.5pi,L_Plot)
    FSS_Plot[1,1] = ax_scal = Axis(fig,xlabel = L"μ",ylabel = L"\xi/L")

    mu_show = (0.,0.3,0.8)

    SqsGFMC = [SW.expand_Sq.(getSq(res_36,tau=15,mu=mu,L=L_Plot)) for mu in mu_show]

    SqMat = mean.(SqsGFMC)
    SqErr = std.(SqsGFMC)
    

    with_theme(theme_PiTicks()) do 
        Sq_Heatmaps[1,1] = ax_mu1 = Axis(fig;xlabel = L"q_x",ylabel = L"q_y",aspect=1,PiTicksArgs...)
        Sq_Heatmaps[1,3] = ax_mu2 = Axis(fig;xlabel = L"q_x",aspect=1,yticklabelsvisible = false,PiTicksArgs...)
        Sq_Heatmaps[1,5] = ax_mu3 = Axis(fig;xlabel = L"q_x",aspect=1,yticklabelsvisible = false,PiTicksArgs...)

    
        linkyaxes!(ax_mu1,ax_mu2,ax_mu3)


        colorrange = extrema(stack(SqMat))

        for (i,ax) in enumerate((ax_mu1,ax_mu2,ax_mu3))
            SqCont = SW.getSqCont(SqMat[i])
            hm = heatmap!(ax,qx,qy,SqCont.(Iterators.product(qx,qy)),colormap = :viridis;
            # colorrange
            )
            Colorbar(Sq_Heatmaps[1,2i],hm,ticks = SimpleTicks(),width = Relative(0.05),height = Relative(0.8))
        end

        # rowsize!(Sq_Heatmaps,1,Relative(0.01))
        # rowsize!(Sq_Heatmaps,2,Relative(0.9999))

    end
    kpath = ["Γ","X","X'","Γ"]
    pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],500)
    kpointlabels = Makie.latexstring.(kpath)
    tRange = eachindex(p1)

    xygrid = [(x,y) for x in qx, y in qy]

    
    
    with_theme(theme_SimpleTicks()) do 
        SqCuts[1,1] = axPath1 = Axis(fig,ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
        )
        SqCuts[1,2] = axPath2 = Axis(fig,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
        )
        SqCuts[1,3] = axPath3 = Axis(fig,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
        )
        # linkyaxes!(axPath1,axPath2,axPath3)

        tRange_new,p1_discrete = rasterCurve(p1,xygrid,tRange)
        for (i,ax,SqMat,SqErr) in zip(eachindex(SqsGFMC), (axPath1, axPath2, axPath3), SqMat, SqErr)
            SqFunc = SW.getSqCont(SqMat)
            SqErrFunc = SW.getSqCont(SqErr)
            Sqcut = [SqFunc(x,y) for (x,y) in xygrid[p1_discrete]]
            Sqerrcut = [SqErrFunc(x,y) for (x,y) in xygrid[p1_discrete]]
            scatter!(ax, tRange_new, Sqcut, marker = '∘', markersize = 10, color = :black)
            errorbars!(ax, tRange_new, Sqcut, Sqerrcut, whiskerwidth = 6, linewidth = 0.5, color = :black)
            
            # Add lines corresponding to the best field theory fit
            fittingCoefs = optimizeCoeffs(SqMat)
            SqFT = [SqFieldTheory(q, fittingCoefs) for q in xygrid[p1_discrete]]
            lines!(ax, tRange_new, SqFT, color = :red, linestyle = :dash)
        end
        SqRandomCuts[1,1] = axPath_Random = Axis(fig,ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),yticks = SimpleTicks()
        )


        text!(axPath_Random,Point(100,1),text="TODO!",color = :black,align = (:center,:center),fontsize = 40)

    end

    with_theme(theme_PiTicks()) do 
        BCorrPlot[1,1] = ax_BCorr = Axis(fig;xlabel = L"q_x",ylabel = L"q_y",aspect=1,PiTicksArgs...)

        Bq = [cos(qx-qy)^2 + cos(qx+qy)^2 for qx in qx, qy in qy]
        hm_B = heatmap!(ax_BCorr,qx,qy,Bq,colormap = Makie.cgrad(:thermal,rev=false))

        text!(ax_BCorr,Point(pi/2,pi/2),text="TODO!",color = :black,align = (:center,:center),fontsize = 40)
    end

    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
        

    # return fig
    # ax2 = insetAtPoint(fig,ax,(0.6,3.2),(110,60))
    # ax2 = Axis(fig[2,1],xlabel = L"μ",ylabel = L"\xi/L")
    # linkaxes!(ax,ax2)

    Linestyles = [:dash,:dot,:solid]
    # allmus = [musSmall,musMedium,mus]
    # allxis = [xisSmall,xisMedium,xis]
    scatterkwargs = Dict(16 => (;marker = '+'),20 => (;marker = '▲'),24 => (;marker = '●' ) ,28 => (;marker = '×',markersize =18),36 => (;marker = '■',markersize =10))
    for (L,linestyle) in zip((20,24,28),Linestyles)
    # for (L,linestyle) in zip(keys(res),Linestyles)
        # Sq = res[L].Sqs[:,:,5,:,:]
        mus = xi_res.mus[L]
        muFilter = findall(x->x<=1,mus)

        xiLs = xi_res.xis[L][muFilter]
        xiLserr = xi_res.xis_err[L][muFilter]
        musPlot = mus[muFilter]
        
        scatterlines!(ax_scal,musPlot,xiLs,label = L"L=%$L";linestyle,scatterkwargs[L]...)
        errorbars!(ax_scal,musPlot,xiLs,xiLserr,whiskerwidth = 5)
    end

    let L=36,linestyle = :dashdot
        musPlot = xi_res.mus[L]
        xiLs = xi_res.xis[L]
        xiLserr = xi_res.xis_err[L]

        scatterlines!(ax_scal,musPlot,xiLs,label = L"L=%$L";linestyle,scatterkwargs[L]...)
        errorbars!(ax_scal,musPlot,xiLs,xiLserr,whiskerwidth = 5)

    end

    inset = insetAtPoint(fig,ax_scal,(0.7,3.8),0.8 .*(110,60),xlabel = L"L",ylabel = L"μ_c")
    scatterlines!(inset,crossings.L,crossings.crossings)

    ylims!(ax_scal,0,5)
    # vlines!(ax_scal,[mu_c],color = :grey,linestyle = :dash)
    # text!(ax_scal,Point(0.22,3.5),text=L"μ_c = %$mu_c",color = :grey,align = (:left,:center))
    # ylims!(ax2,0.04,1.8)
    axislegend(ax_scal,position = :lb)

    rowsize!(fig.layout,1,Relative(0.3))
    rowsize!(fig.layout,2,Relative(0.3))
    colgap!(Sq_Heatmaps,1,-40)
    colgap!(Sq_Heatmaps,2,-40)
    colgap!(Sq_Heatmaps,3,-40)
    colgap!(Sq_Heatmaps,4,-40)
    colgap!(Sq_Heatmaps,5,-40)

    fig
    
end

##

with_theme(theme_SimpleTicks()) do 
    L = 36
    muIndex = findfirst(>=(0.8),res_36.mu)

    μ = res_36.mu[muIndex]

    SqsGFMC = SW.expand_Sq.(getSq(res_36,tau=26,mu=μ,L=L))
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
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    
    SqFT = [SqFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]
    heatmap!(axFT,qx,qy,SqFT)
    q_path(r,phi) = (r*cos(phi),r*sin(phi))
    qr = LinRange(0,.35pi,100)
    
    colors = (:red,:blue,:magenta)
    

    colorFT = :black
    colorGFMC = :red

    kpath = ["Γ","X","X'","Γ"]
    pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],500)
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

    # for (phi,color) in zip([0,pi/4],colors)
    #     qpoints_raw = q_path.(qr,phi)
    #     qpoints = sort!(unique!(roundToTrueMomenta.(qpoints_raw,size(SqMat,1)-1)), by = SW.norm)

    #     Sqcut = Sq.(qpoints)
    #     Sqerrcut = Sqerr.(qpoints)
        
    #     # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
    #     SqFT = [SqFieldTheory(q,fittingCoefs...) for q in qpoints]
    #     scatter!(ax,qpoints,marker = '×' ,color = color)
    #     scatterlines!(axFT,Point.(qpoints),color = color,linestyle = :dash,marker = '●',markersize = 4)
    #     qnorms_sq = SW.norm.(qpoints).^2
    #     scatter!(ax2,qnorms_sq,Sqcut,
    #     marker = '×',markersize = 15,color = color)
    #     errorbars!(ax2,qnorms_sq,Sqcut,Sqerrcut,color = color,whiskerwidth = 6,linewidth=0.5)
    #     scatterlines!(ax2,qnorms_sq,SqFT,color = color,linestyle = :dash,marker = '●',markersize = 4)
    # end
    rowsize!(fig.layout,1,Relative(0.5))
    text!(axPath,Point(tRange[end-end÷8],1.2),text=L"μ = %$μ",align = (:center,:center))
    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    # Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end
