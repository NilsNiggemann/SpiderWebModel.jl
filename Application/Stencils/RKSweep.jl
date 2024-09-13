using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
include("plottingUtils.jl")
##
function prepResults(folder,mufilter)
    files = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
    filter!(contains("mu=$(mufilter)_"),files)
    return files
    res = vcat(SW.readResults.(files,5000,discardStart = 5000)...)
    energies = [h5read(file,"energies") for file in files]
    TotalWeights = [h5read(file,"TotalWeights") for file in files]
    # mus = [h5read(file,"mu") for file in files]
    # taus = [h5read(file,"tau") for file in files]

    ens = stack(SW.getEnergies.(TotalWeights,energies,1,1000))
end
# entest = prepResults("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equiv_open/",0.30)
entest = prepResults("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/L=32/",0.65)
##
files = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=32/") for file in files]#[1:2:end]
# files = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_equil/eval/L=30/") for file in files]#[1:2:end]
# filter!(!contains("mu=-0.18"),files)
# filter!(!contains("mu=-0.2"),files)
##
energies = stack([h5read(file,"energies") for file in files])
mus = [h5read(file,"mu") for file in files]
Sqs = stack([h5read(file,"SqsGFMC/200") for file in files])
taus = [h5read(file,"tau") for file in files]
##
with_theme(theme_SimpleTicks()) do
    ind = 3
    Nsites = length(Sqs[1:end-1,1:end-1,1,ind])
    enmean = mean(energies,dims=2)[1:250,1,ind] ./ Nsites

    enstd = std(energies,dims=2)[1:250,1,ind] ./ Nsites
    # enmean = mean(entest,dims=2)[:,1]
    # enstd = std(entest,dims=2)[:,1]
    tau = taus[ind]
    lines((eachindex(enmean).-1) .*tau,enmean,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    band!((eachindex(enmean).-1) .*tau,enmean - enstd , enmean + enstd,color = (:black,0.2))
    current_figure()

end
##
with_theme(theme_SimpleTicks()) do 
    tauindices = [round(Int,20 ÷ tau) for tau in taus]
    energies_slice = zeros(size(energies,2),size(energies,3))
    for (i,tau) in enumerate(tauindices)
        energies_slice[:,i] .= @view energies[tau,:,i]
    end

    enmean = mean(energies_slice,dims=1)[1,:]
    enstd = std(energies_slice,dims=1)[1,:]
    Nsites = 30^2
    push!(enmean,0)
    push!(enstd,0)
    mus2 = copy(mus)
    push!(mus2,1)
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
    fig
end

##
with_theme(theme_PiTicks()) do 
    Sq = mean(Sqs,dims=3)[:,:,1,:] ./ 4
    # muPlot = [-0.06,0.2,0.3,0.6,0.94,1.1]
    muPlot = [0.4,0.6,0.9,1.05]

    fig = Figure(fontsize = 22,size = 200 .*(length(muPlot),1.4))
    ticks = PiTicks([0,pi])

    # ax1 = Axis(fig[1,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = %$(mus[3])")
    # SqMat = Sq[:,:,3]
    # SqFunc = SW.getSqCont(SqMat)

    kx = ky = trueMomenta(-pi/2,1.5pi,size(Sq,1)-1)
    # Sqpl = SqFunc.(Iterators.product(kx,ky))
    # hm = heatmap!(ax1,kx,ky,Sqpl,colormap = :viridis)
    # muPlot = [0.9,0.92,0.94,0.96]
    mupls = mus[[findfirst(>=(mu),mus) for mu in muPlot]]
    axes = [Axis(fig[1,i],aspect=1,title = L"μ = %$(mupls[i])",yticklabelsvisible=i==1,xticks=ticks,yticks=ticks,xlabel = L"q_x",ylabel = L"q_y", ylabelvisible = i==1) for i in eachindex(muPlot)]

    for (i,ax) in enumerate(axes)
        i_mu = findfirst(>(muPlot[i]),mus)
        SqMat = Sq[:,:,i_mu]
        mupl = mus[i_mu]
        SqFunc = SW.getSqCont(SqMat)
        Sqpl = SqFunc.(Iterators.product(kx,ky))
        heatmap!(ax,kx,ky,Sqpl,colormap = :viridis)
    end
    fig
end
##
muIndex = 10
SqsGFMC = Sqs[:,:,:,muIndex]./ 4
SqMat = dropdims(mean(SqsGFMC,dims=3),dims=3)
SqErr = dropdims(std(SqsGFMC,dims=3),dims=3)
fittingCoefs = optimizeCoeffs(SqMat)
##
with_theme(theme_SimpleTicks()) do 

    μ = mus[muIndex]
    fig = Figure(size = 120 .* (4,5))

    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks(), yticks = PiTicks())

    ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks(), yticks = PiTicks(),ylabelvisible = false,yticklabelsvisible = false)

    ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
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
    
        
        axPath = Axis(fig[3,1:2],ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
        )
        tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
        
    
        p1_points = xygrid[p1_discrete]

        Sqcut = [Sq(x,y) for (x,y) in p1_points]
        Sqerrcut = [Sqerr(x,y) for (x,y) in p1_points]
        SqFT = [SqFieldTheory(q,fittingCoefs...) for q in p1_points]
        
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

    for (phi,color) in zip([0,pi/4],colors)
        qpoints_raw = q_path.(qr,phi)
        qpoints = sort!(unique!(roundToTrueMomenta.(qpoints_raw,size(SqMat,1)-1)), by = SW.norm)

        Sqcut = Sq.(qpoints)
        Sqerrcut = Sqerr.(qpoints)
        
        # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
        SqFT = [SqFieldTheory(q,fittingCoefs...) for q in qpoints]
        scatter!(ax,qpoints,marker = '×' ,color = color)
        scatterlines!(axFT,Point.(qpoints),color = color,linestyle = :dash,marker = '●',markersize = 4)
        qnorms_sq = SW.norm.(qpoints).^2
        scatter!(ax2,qnorms_sq,Sqcut,
        marker = '×',markersize = 15,color = color)
        errorbars!(ax2,qnorms_sq,Sqcut,Sqerrcut,color = color,whiskerwidth = 6,linewidth=0.5)
        scatterlines!(ax2,qnorms_sq,SqFT,color = color,linestyle = :dash,marker = '●',markersize = 4)
    end
    rowsize!(fig.layout,1,Relative(0.4))

    Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end

##
res1 = SW.readResults(entest[2],2000)[1:2]
Sqs = SW.getSqsGFMC(res1,200;discardborder = 4)
##
with_theme(theme_PiTicks()) do 
    Sq = mean(Sqs) ./ 4

    fig = Figure(fontsize = 22,size = 400 .*(1.4,1.4))
    ticks = PiTicks([0,pi])

    ax1 = Axis(fig[1,1],aspect = 1,xlabel = L"q_x",ylabel = L"q_y",xticks = ticks,yticks = ticks,title = L"μ = 0.6")
    kx = ky = trueMomenta(-pi/2,1.5pi,size(Sq,1)-1)
    SqFunc = SW.getSqCont(Sq)
    heatmap!(ax1,kx,ky,SqFunc.(Iterators.product(kx,ky)),colormap = :viridis)
    fig
end
##
mags = SW.getObs(res1[1],float,100)

##
let
    mags1 = mags[begin+4:end-4,begin+4:end-4] ./2
    # mags1 = mags./2
    # for I in CartesianIndices(mags1)
    #     i,j = Tuple(I)
    #     if iseven(i+j)
    #         mags1[I] = NaN
    #     end
    # end
    fig,ax,hm = heatmap(mags1)
    Colorbar(fig[1,2], hm)
    fig
end