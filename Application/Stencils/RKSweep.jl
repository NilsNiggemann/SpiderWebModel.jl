import SpiderWebModel as SW
using CairoMakie, MakieHelpers,Statistics
using SpiderWebModel.HDF5
using DataFrames
include("plottingUtils.jl")

##
function getxi(Sq,I::CartesianIndex,dI::CartesianIndex)
    L = size(Sq,1)-1
    if Sq[I] < Sq[I+dI]
        return 0
    end
    xi_L = sqrt(Sq[I]/Sq[I+dI] -1 )
    return xi_L * L
end

function getxi(Sq,I::CartesianIndex)
    # neighbors = [(-1,-1)]
    neighbors = [(i,j) for i in -1:1,j in -1:1 if i != 0 || j != 0]
    # neighbors = [(1,0),(0,1),(-1,0),(0,-1)]
    return maximum(getxi(Sq,I,CartesianIndex(dI)) for dI in neighbors)
end

function getxi(Sq)
    # L = size(Sq,1)-1
    I = argmax(Sq)
    # return Sq[I] / L^2
    return getxi(Sq,I)
end

##
function getXis(Sqs)
    xis = zeros(size(Sqs,4))
    for (i,ii) in enumerate(axes(Sqs,4))
        Sq = @views dropmean(Sqs[:,:,:,ii],dims=3)

        xis[i] = getxi(Sq./4)
    end
    return xis
end
function getXis_err(Sqs,I)
    xis = zeros(size(Sqs)[3:4])
    for (i,ii) in enumerate(axes(Sqs,4))
        for (j,jj) in enumerate(axes(Sqs,3))
            Sq = @views Sqs[:,:,jj,ii]

            xis[j,i] = getxi(Sq./4,I)
        end
    end
    ximean = dropmean(xis,dims=1)
    xistd = dropstd(xis,dims=1)
    return (;ximean,xistd)
end

function getXis(Sqs::AbstractVector{<:AbstractMatrix},I)
    xis = zeros(length(Sqs))
    for (i,Sq) in enumerate(Sqs)
        xis[i] = getxi(Sq./4,I)
    end
    return xis
end

##
function getFiles(L)
    files = [joinpath(root,file) for (root,_,files) in walkdir( "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=$L/") for file in files]
end
function prepResults(folder,mufilter)
    files = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
    filter!(contains("mu=$(mufilter)_"),files)
    # return files
    res = vcat(SW.readResults.(files,5000)...)
    return res
    energies = [h5read(file,"energies") for file in files]
    TotalWeights = [h5read(file,"TotalWeights") for file in files]
    # mus = [h5read(file,"mu") for file in files]
    # taus = [h5read(file,"tau") for file in files]

    ens = stack(SW.getEnergies.(TotalWeights,energies,1,1000))
end
# entest = prepResults("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equiv_open/",0.30)
# entest = prepResults("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/L=32/",0.65)
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

function getRes(folder)
    files = let
        filesunsrt = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
        fileformats = file_format.(filesunsrt)
        invalid_files = findall(iszero,fileformats)
        if !isempty(invalid_files)
            println("invalid files:")
            println(filesunsrt[invalid_files])
        end
        filesunsrt = [f for (f,i) in zip(filesunsrt,fileformats) if i == 1]
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
    Sqs = stack([h5read(file,"SqsGFMC") for file in files])
    taus = [h5read(file,"tau") for file in files]
    p_Sq = stack([h5read(file,"p_Sq") for file in files])
    return (;energies,mus,Sqs,taus,files,p_Sq)
end

function getSq_tau(res,tau)
    p_Sq = res.p_Sq
    Dtaus = reshape(res.taus,1,size(p_Sq,2))
    taus = p_Sq .* Dtaus
    tauInds = [findfirst(>=(tau),t) for t in eachcol(taus)]
    # return tauInds
    return stack(res.Sqs[:,:,ti,:,i] for (i,ti) in enumerate(tauInds))
end

function getSq_tau(res::DataFrame,tau)
    Dtaus = res.tau

    taus = [axes(Sqs,3) .* t for (Sqs,t) in zip(res.Sq,Dtaus)]
    tauInds = [findfirst(>=(tau),t) for t in taus]

    # return tauInds
    return [res.Sq[i][:,:,ti] for (i,ti) in enumerate(tauInds)]
end
function getSq_tau(res::DataFrame,::Nothing)
    return res.Sq
end

function getSq(res::DataFrame;tau,mu,L)
    res_new = filter(row -> row.mu == mu && row.L == L,res)
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

res = Dict(L=>getRes("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=$L/") for L in (16,20,24,28))

##

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
res_36 = getRes_2("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=36/")
##
with_theme(theme_SimpleTicks()) do
    ind = 1
    L = 28
    Nsites = length(res[28].Sqs[1:end-1,1:end-1,1,ind,begin])
    enmean = mean(res[28].energies,dims=2)[1:250,1,ind] ./ Nsites

    enstd = std(res[28].energies,dims=2)[1:250,1,ind] ./ Nsites
    # enmean = mean(entest,dims=2)[:,1]
    # enstd = std(entest,dims=2)[:,1]
    tau = res[28].taus[ind]
    errlines((eachindex(enmean).-1) .*tau,enmean,enstd,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # lines((eachindex(enmean).-1) .*tau,enmean,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # band!((eachindex(enmean).-1) .*tau,enmean - enstd , enmean + enstd,color = (:black,0.2))
    # current_figure()

end
##
with_theme(theme_SimpleTicks()) do 
    L = 28

    tauindices = [round(Int,8 ÷ tau) for tau in res[L].taus]
    energies_slice = zeros(size(res[L].energies,2),size(res[L].energies,3))
    for (i,tau) in enumerate(tauindices)
        energies_slice[:,i] .= @view res[L].energies[tau,:,i]
    end

    enmean = dropdims(mean(energies_slice,dims=1),dims=1)
    enstd = dropdims(std(energies_slice,dims=1),dims=1)
    Nsites = L^2
    # push!(enmean,0)
    # push!(enstd,0)
    mus2 = copy(res[L].mus)
    # push!(mus2,1)
    ord = sortperm(mus2)
    mus2 = mus2[ord]
    enmean = enmean[ord] ./ Nsites
    enstd = enstd[ord] ./ Nsites

    fig = Figure(fontsize = 22,size = (800,400))

    axDE = Axis(fig[1, 1], xlabel = L"μ", ylabel = L"dE/d\mu/N_\text{sites}",
        yaxisposition=:right,yticklabelcolor=:red,
        yticks = SimpleTicks(),
        xlabelvisible = false,
        ygridvisible=false,
        xgridvisible=false,
        xticklabelsvisible= false,
        xticksvisible= false,
        # xminorticksvisible=false
        ylabelvisible=true,
        ylabelcolor = :red,
        yminorticksvisible=true,yminorticks = IntervalsBetween(4),
    )
    ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"E/N_\text{sites}",
    xminorticksvisible=true,xticks = -0.2:0.2:1.2,xminorticks = IntervalsBetween(2),yminorticksvisible=true,yminorticks = IntervalsBetween(4),
    )
    plaqPhase = [-1,0.16]
    SLphase = [0.16,1]
    OrderedPhase = [1,3]
    band!(ax,plaqPhase,[-300,-300],[10,10],color = (:blue,0.2))

    band!(ax,SLphase,[-300,-300],[10,10],color = (:green,0.2))
    band!(ax,OrderedPhase,[-300,-300],[10,10],color = (:red,0.2))

    # band!(axDE,,dEdmu .- std(dEdmu),dEdmu .+ std(dEdmu),color = (:red,0.2))


    linkxaxes!(ax,axDE)


    dEdmu = diff(enmean)./diff(mus2)

    # dEdmus = diff(energies,dims=3) ./ diff(mus)
    # enps = eachrow(energies[1000,:,:])
    # return enps
    # dE = diff(enp,dims=2)
    # dEdmus = stack(eachrow(dE) ./ diff(mus))
    # dEdmu = mean(dEdmus,dims=2)[:]
    # return dEdmu
    scatterlines!(axDE,mus2[2:end],dEdmu,label = "dE/dμ",color = :red,marker = :rect,linestyle = :dash)
    scatterlines!(ax,mus2,enmean,label = "Energy",color = :black)
    errorbars!(ax,mus2,enmean,enstd,label = "std error",color = :black,whiskerwidth = 8)
    vlines!(ax,[0.16,1],color = :grey,linestyle = :dash)
    xlims!(ax,extrema(mus2)...)
    ylims!(ax,extrema(enmean)...)
    # return enmean, mus2
    ylims!(axDE,extrema(dEdmu)...)
    fig
end

##
with_theme(theme_PiTicks()) do
    L = 28
    Sq = dropmean(res[L].Sqs,dims=4)[:,:,10,:] ./ 4
    # muPlot = [-0.06,0.2,0.3,0.6,0.94,1.1]
    # muPlot = [0.0,0.4,0.9,1.05]
    muPlot = [0.2,0.6,]

    fig = Figure(fontsize = 22,size = 200 .*(length(muPlot),1.4))
    ticks = PiTicks([0,pi])

    # ax1 = Axis(fig[1,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = %$(mus[3])")
    # SqMat = Sq[:,:,3]
    # SqFunc = SW.getSqCont(SqMat)

    kx = ky = trueMomenta(-pi/2,1.5pi,size(Sq,1)-1)
    # Sqpl = SqFunc.(Iterators.product(kx,ky))
    # hm = heatmap!(ax1,kx,ky,Sqpl,colormap = :viridis)
    # muPlot = [0.9,0.92,0.94,0.96]
    mupls = res[L].mus[[findfirst(>=(mu),res[L].mus) for mu in muPlot]]
    spinconf = SW.SpinConfig(SW.periodicStateLoops(8),1)
    ax0 = Axis(fig[1,0];SW.getConfigAxis(spinconf)...,xticks = 1:2:8 ,yticks = 1:2:8,xlabel = L"x",ylabel = L"y")
    
    SW.plotSpinConfig!(ax0,spinconf)

    axes = [Axis(fig[1,i],aspect=1,title = L"μ = %$(mupls[i])",yticklabelsvisible=i==1,xticks=ticks,yticks=ticks,xlabel = L"q_x",ylabel = L"q_y", ylabelvisible = i==1) for i in eachindex(muPlot)]

    for (i,ax) in enumerate(axes)
        i_mu = findfirst(>(muPlot[i]),res[L].mus)
        SqMat = Sq[:,:,i_mu]
        mupl = res[L].mus[i_mu]
        SqFunc = SW.getSqCont(SqMat)
        Sqpl = SqFunc.(Iterators.product(kx,ky))
        heatmap!(ax,kx,ky,Sqpl,colormap = :viridis)
    end
    colgap!(fig.layout,2,0)
    # colgap!(fig.layout,3,0)
    # colgap!(fig.layout,4,0)
    Label(fig[1,0, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    Label(fig[1,1, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
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

with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 22,size = 400 .*(1.4,1.))
    ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"\xi/L")
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
        Sq = getSq_tau(res[L],10)
        k = trueMomenta(0,2pi,size(Sq,1)-1)
        i_k = findfirst(==(pi/2),k)
        xis = getXis_err(Sq,CartesianIndex(i_k,i_k))
        mus = res[L].mus
        muFilter = findall(x->x<=1,mus)

        xiLs = xis.ximean[muFilter] ./ L
        xiLserr = xis.xistd[muFilter] ./ L
        musPlot = mus[muFilter]
        scatterlines!(ax,musPlot,xiLs,label = L"L=%$L";linestyle,scatterkwargs[L]...)
        errorbars!(ax,musPlot,xiLs,xiLserr,whiskerwidth = 5)
    end
    let L=36,linestyle = :dashdot
        unique_mus = unique(res_36.mu)
        k = trueMomenta(0,2pi,L)
        i_k = findfirst(==(pi/2),k)

        Sq = [getSq(res_36,mu = mu,L = L,tau=10) for mu in unique_mus]

        xis = [getXis(Sq,CartesianIndex(i_k,i_k)) for Sq in Sq]

        xiLs = [mean(xi) for xi in xis] ./ L
        xiLserr = [std(xi) for xi in xis] ./ L

        musPlot = unique_mus
        scatterlines!(ax,musPlot,xiLs,label = L"L=%$L";linestyle,scatterkwargs[L]...)
        errorbars!(ax,musPlot,xiLs,xiLserr,whiskerwidth = 5)

    end


    ylims!(ax,0,5)
    mu_c = 0.18
    vlines!(ax,[mu_c],color = :grey,linestyle = :dash)
    text!(ax,Point(0.22,3.5),text=L"μ_c = %$mu_c",color = :grey,align = (:left,:center))
    # ylims!(ax2,0.04,1.8)
    axislegend(ax,position = :lb)
    fig
    
end
##
a = let L=36
    unique_mus = unique(res_36.mu)
    # xis = [getxi(mean(filter(x->x.mu == mu,res_36).Sq)) for mu in unique_mus]
    Sqs = [getSq(res_36,mu = mu,L = L,tau=10) for mu in unique_mus]
end
##
KPoints = Dict([
    "Γ" => SVector(0,0),
    "X" => SVector(pi,0),
    "M" => SVector(pi,pi),
    "X'" => SVector(0,pi)
    ])



with_theme(theme_SimpleTicks()) do 
    L = 28
    muIndex = findfirst(>=(0.8),res[L].mus)
    SqsGFMC = getSq_tau(res[L],6)[:,:,:,muIndex]./ 4
    SqMat = dropmean(SqsGFMC,dims=3)
    SqErr = dropstd(SqsGFMC,dims=3)
    fittingCoefs = optimizeCoeffs(SqMat)
    
    μ = res[L].mus[muIndex]
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
