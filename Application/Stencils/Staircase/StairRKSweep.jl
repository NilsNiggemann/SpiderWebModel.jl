import SpiderWebModel as SW
using CairoMakie, MakieHelpers,Statistics
using SpiderWebModel.HDF5
using DataFrames
using MakieExtra
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

function getEnergy_tau(res::DataFrame,tau;kwargs...)
    res = get_res(res;kwargs...)
    Dtaus = res.tau

    taus = [axes(Energies,1) .* t for (Energies,t) in zip(res.Energy,Dtaus)]
    tauInds = [findfirst(>=(tau),t) for t in taus]

    # return tauInds
    return [res.Energy[i][ti] for (i,ti) in enumerate(tauInds)]
end

function getRes_2(folder)
    files = let
        filesunsrt = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
        fileformats = file_format.(filesunsrt)
        invalid_files = findall(iszero,fileformats)
        if !isempty(invalid_files)
            # rm.(filesunsrt[invalid_files])
            # println("invalid files:")
            # println(filesunsrt[invalid_files])
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
    NWalkers = [h5read(file,"NWalkers") for file in files]

    res= DataFrame(;L,mu,Energy,Sq,tau,NWalkers,files)
    sort!(res,[:mu,:L])
end
##
res = getRes_2(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/StairCase_merge/")
# res2 = getRes_2(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/StairCase_2/")
# res3 = getRes_2(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/StairCase_3/")
# res2fine = getRes_2(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/StairCase_fine_2/")
# res2fine_2 = getRes_2(ENV["MYSCRATCH"]*"Spiderweb/DataS1_CT_RK_equil/StairCase_fine/")
# filter!(x->!(x.mu in res2.mu && x.L in res2.L),res)
# filter!(x->!(x.mu in res3.mu && x.L in res3.L),res)
# res = vcat(res,res2,res2fine,res3)
# filter!(x->!(x.mu ==0.8 && x.L==36 && x.NWalkers <20000),res)
# filter!(x->x.mu>=0.6,res)
sort!(res,[:mu,:L])

##
with_theme(theme_SimpleTicks()) do
    L = 36
    resL = get_res(res,L=L,mu=0.8)

    Nsites = length(getSq_tau(res,8.))
    enmean = mean(resL.Energy)./ Nsites

    enstd = std(resL.Energy)./ Nsites
    # enmean = mean(entest,dims=2)[:,1]
    # enstd = std(entest,dims=2)[:,1]
    tau = only(unique(resL.tau)) .* (eachindex(enmean) .-1)
    errlines(tau,enmean,enstd,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # current_figure()

end
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22, size = 200 .* (3, 2))
    # ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"E/N_\textrm{sites}", xticks = SimpleTicks(), yticks = SimpleTicks())
    ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"E/N_\textrm{sites} /(1-\mu)", xticks = SimpleTicks(), yticks = SimpleTicks())

    Linestyles = [:dash, :dot, :dashdot, :solid, :solid]
    scatterkwargs = Dict(
        28 => (;marker = '●'),
        32 => (;marker = '▼', markersize = 10),
        36 => (;marker = '▲', markersize = 10),
    )
    inset_En = Axis(fig[1,1],width = Relative(0.3),height = Relative(0.5),halign=0.2, valign=0.85)

    for (L, linestyle) in zip((28,32,36), Linestyles)
        resL = get_res(res, L=L)
        filter!(x->x.mu<=0.9,resL)
        mus = unique(resL.mu)
        energies = [getEnergy_tau(resL,10;mu) for mu in mus]

        enmean = mean.(energies) ./ L^2  ./ (1 .-mus)
        enstd = std.(energies) ./ L^2  ./ (1 .-mus)

        scatterlines!(ax, mus, enmean, label=L"L=%$L"; linestyle, scatterkwargs[L]...)
        errorbars!(ax, mus, enmean, enstd, whiskerwidth=5)

        mufilter = findall(x -> 0.85 >= x>= 0.78 ,mus)
        scatterlines!(inset_En, mus[mufilter], enmean[mufilter], label=L"L=%$L"; linestyle, scatterkwargs[L]...)
        errorbars!(inset_En, mus[mufilter], enmean[mufilter], enstd[mufilter], whiskerwidth=5)
    end

    # inset_En2 = Axis(fig[1,1],width = Relative(0.3),height = Relative(0.3),halign=0.8, valign=0.2)

    # translate!(inset_En.blockscene,0,0,1000)
    # translate!(inset_En2.blockscene,0,0,1000)
    zoom_lines!(ax,inset_En)

    axislegend(ax, position=:rb)
    fig
end
##
with_theme(theme_SimpleTicks()) do
    L = 36
    resL = get_res(res,mu=0.81,L=L)

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
        range = 1:120
        # scatterlines!(ax,tau[range],SqMat[i,j,range],marker = '×')
        # errorbars!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],whiskerwidth = 6,linewidth=0.5)
        errlines!(ax,tau[range],SqMat[i,j,range],SqErr[i,j,range],linewidth=0.5)
    end
    fig
end

##
xis = getXiCol(res,(pi/2,pi/2);tau = 3)
Sqmax = getSqMaxCol(res,(pi/2,pi/2),tau = 3)
res.xi = xis
res.Sqmax = Sqmax

with_theme(theme_PiTicks()) do
    LPlot = 36
    # muPlot = [0.8,0.9,0.9,1.0]
    muPlot = [0.8,0.81,0.9,1.0]
    
    fig = Figure(fontsize = 22,size = 300 .*(2,2.5))

    fig_xi = GridLayout()
    fig.layout[1,1:2] = fig_xi
    
    # axis = Axis(fig_xi[1,1],xlabel = L"μ",ylabel = L"\xi/L",xticks = SimpleTicks(),yticks = SimpleTicks())
    axis = Axis(fig_xi[1,1],xlabel = L"μ",ylabel = L"\mathcal{S}(\frac{\pi}{2},\frac{\pi}{2})",xticks = SimpleTicks(),yticks = SimpleTicks())

    for L in filter!(<(40), unique(res.L))
        resFilt = filter(x->x.L==L,res)
        unique_mus = unique(resFilt.mu)
        data = map(unique_mus) do mu
            resmu = filter(x->x.mu==mu,resFilt)
            # return mean(resmu.xi)/L, std(resmu.xi)/L ,mu
            # return mean(resmu.Sqmax), std(resmu.Sqmax) ,mu
            return mean(resmu.Sqmax), std(resmu.Sqmax) ,mu
        end
        xi = getindex.(data,1)
        xierr = getindex.(data,2)
        mu = getindex.(data,3)
        scatterlines!(axis,mu,xi,marker = '×',markersize = 10,label = "L = $L")
        errorbars!(axis,mu,xi,xierr,whiskerwidth = 6,linewidth=0.5)
    end
    # xlims!(axis,0.6,1)
    # ylims!(axis,0.0,70)
    axislegend(axis,position = :rt)
    fig_SQ = GridLayout()
    fig.layout[2,1:2] = fig_SQ


    ticks = PiTicks([0, pi])

    axInds = [(i,(2j-1)) for (i, j) in Iterators.product(1:2, 1:2)]
    ColbarInds = [(i,j+1) for (i, j) in axInds]

    axes = [Axis(fig_SQ[i, j], aspect=1,
    yticklabelsvisible=j==1, xticks=ticks, yticks=ticks, xlabel=L"q_x", ylabel=L"q_y", ylabelvisible=j==1,xticklabelsvisible= i==2, xlabelvisible = i==2) 
    for (i, j) in axInds]
    

    for (idx, mu) in enumerate(muPlot)
        i, j = divrem(idx - 1, 2) .+ 1
        SqMat = SW.expand_Sq(mean(getSq(res, tau=12., L=LPlot, mu=mu)))

        kx = ky = trueMomenta(-pi/2, 1.5pi, LPlot)

        SqFunc = SW.getSqCont(SqMat)
        Sqpl = SqFunc.(Iterators.product(kx, ky))
        hm = heatmap!(axes[i, j], kx, ky, Sqpl, colormap=:viridis)
        Colorbar(fig_SQ[ColbarInds[idx]...], hm, ticks=SimpleTicks(), vertical=true, flipaxis=true, label=L"\mathcal{S}(\mathbf{q})", width=Relative(0.5),labelvisible=ColbarInds[idx][2] == 4)
        
        text!(axes[i, j], Point(0.5,0.99), text=L"\mu = %$(mu)", align=(:center, :top),space = :relative, strokecolor = (:black,0.3), strokewidth = 4)
        text!(axes[i, j], Point(0.5,0.99), text=L"\mu = %$(mu)", align=(:center, :top),space = :relative,color = :white)

        # band!(axes[i, j], [-0.5pi, 0.5pi], [1.1pi, 1.1pi], [1.5pi, 1.5pi], color=(:black, 0.5))
        # text!(axes[i, j], Point(0, 1.3pi), text=L"μ = %$mu", color=:white, align=(:center, :center))
    end
    rowsize!(fig.layout, 1, Relative(0.4))
    # rowgap!(fig.layout, 1, 5)
    colsize!(fig_SQ, 2, Relative(0.05))
    colsize!(fig_SQ, 4, Relative(0.05))
    Label(fig[1, 1, TopLeft()], L"(a)$$", padding = (-60, 0, -10, 0))
    Label(fig[2, 1, TopLeft()], L"(b)$$", padding = (-60, 0, -10, 0))
    # colgap!(fig.layout, 2, 0)
    save("../../figs/PaperFigs/StairCaseSpin1Overview.pdf", fig)
    fig
end
##

function errorBarLegend(size = 0.5;linekwargs = (;),markerkwargs = (;),kwargs...)
    center = 0.5
    ymin = center - size/2
    ymax = center + size/2

    errorbar = [Point2f(center, ymin), Point2f(center, ymax)]
    [
        LineElement(linepoints = errorbar;kwargs...,linekwargs...),
        MarkerElement(points = errorbar, marker = :hline, markersize = 10;kwargs...,linekwargs...),
        MarkerElement(points = [Point2f(center, center),], marker = '●', markersize = 7;kwargs...,markerkwargs...)
    ]
end

with_theme(theme_SimpleTicks()) do 
    mu = 0.9
    L = 36
    # SqMat = SW.expand_Sq(mean(getSq(res, tau=12., L=L, mu=mu)))
    SqsGFMC = SW.expand_Sq.(getSq(res,tau=12,mu=mu,L=L))
    SqMat = mean(SqsGFMC)
    SqErr = std(SqsGFMC)
    
    # SqErr = dropstd(SqsGFMC,dims=4)[:,:,end,:]
    fittingCoefs = optimizeCoeffsAsym(SqMat)
    # fittingCoefs = [1,2,3.]
    μ = mu
    fig = Figure(size = 140 .* (4,4),fontsize = 22)

    xticks = yticks = PiTicks()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks)

    axFT = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = false,yticklabelsvisible = false)

    # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
    Sq = SW.getSqCont(SqMat)
    Sqerr = SW.getSqCont(SqErr)
    qx = qy = trueMomenta(-pi,2pi,size(SqMat,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    
    SqFT = [AsymFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]

    heatmap!(axFT,qx,qy,SqFT)
    q_path(r,phi) = (r*cos(phi),r*sin(phi))
    qr = LinRange(0,.35pi,100)
    
    colors = (:red,:blue,:magenta)
    

    colorFT = :red
    colorGFMC = :black

    # kpath = ["Γ","X","X'","Γ"]
    kpath = ["D","X","X'","Γ","D"]

    KPointsNew = Dict([
        "Γ" => SVector(0,0),
        "X" => SVector(pi,0),
        "M" => SVector(pi,pi),
        "X'" => SVector(0,pi),
        "D" => SVector(0,-pi),
        ])

    pointlabels,p1 = fetchKPath([KPointsNew[k] for k in kpath],100)
    kpointlabels = Makie.latexstring.(kpath)
    tRange = eachindex(p1)
    xygrid = [(x,y) for x in qx, y in qy]

    
    axPath = Axis(fig[2,1:2],ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),yticks = [0,0.5,1,1.5],
    )
    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    

    p1_points = xygrid[p1_discrete]

    Sqcut = [Sq(x,y) for (x,y) in p1_points]
    Sqerrcut = [Sqerr(x,y) for (x,y) in p1_points]
    SqFT = [AsymFieldTheory(q,fittingCoefs) for q in p1_points]
    # SqFT = [AsymFieldTheory(q,1,10) for q in qpoints]
    lines!(axFT,p1_points,color = colorFT,linestyle = :dash,linewidth = 3)
    # scatterlines!(axFT,p1_points,color = colorFT,linestyle = :dash,marker = '●',markersize = 2)
    # tRange = SW.norm.(p1).^2
    lines!(axPath,tRange,SqFT,color = colorFT,linestyle = :dash,label = L"field theory$$",linewidth = 2)
    for point in kpath
        P = KPointsNew[point]

        align = (:center,:center)
        if point == "D"
            align = (:center,:bottom)
        end
        text!(axFT,Point(P);text=Makie.latexstring(point),align,strokecolor=(:black,0.3),strokewidth=4)
        text!(axFT,Point(P);text=Makie.latexstring(point),align,color = :white,)
    end
    A_fit, r_fit, p_fit = strd.(fittingCoefs )

    # text!(axPath,Point(0.7,0.85),text = L"A = %$(A_fit),\ r = %$(r_fit)\ p = %$(p_fit)",align = (:bottom,:right),space = :relative,color = :black,fontsize = 16)
    markersizeGFMC = 23
    scatter!(axPath,tRange,Sqcut,
    marker = '∘',markersize = markersizeGFMC,color = colorGFMC,label = L"GFMC$$")
    errorbars!(axPath,tRange,Sqcut,Sqerrcut,color = colorGFMC,whiskerwidth = 6,linewidth=0.5,label = L"GFMC$$")
    
    Legend(fig[2, 2], [
        [
        errorBarLegend(0.6,;color = colorGFMC,markerkwargs = (;marker='∘',markersize = markersizeGFMC)),
    ],  
    [LineElement(linepoints = [Point2f(0., 0.5), Point2f(1, 0.5)],color = colorFT,linestyle = :dash,linewidth = 2)
    ]], [L"GFMC$$",L"$A = %$(A_fit),\ r = %$(r_fit)\ p = %$(p_fit)$"], tellheight = false, tellwidth = false,halign = :right, valign = :top,margin = (0,0,0,0),labelsize = 20,
    )
    ylims!(axPath,0,1.75)
    # axislegend(axPath,position = :lt,merge = true)
    rowsize!(fig.layout,1,Relative(0.5))
    rowgap!(fig.layout,1,-10)
    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    Label(fig[1,1, TopLeft()],L"(a)$$",padding = (-60,0,-20,0))
    Label(fig[1,2, TopLeft()],L"(b)$$",padding = (-20,0,-20,0))
    Label(fig[2,1, TopLeft()],L"(c)$$",padding = (-60,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))
    save("../../figs/PaperFigs/StairCaseSpin1_Sq_FT.pdf", fig)
    fig
end

