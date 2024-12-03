using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
include("../plottingUtils.jl")
strd(x;kwargs...) = string(round(x,digits = 1);kwargs...)

function is_valid_file(filename)
    allkeys = ["energies","mu","tau","SqsGFMC","p_Sq"]
    h5open(filename,"r") do file
        all(k->haskey(file,k),allkeys) || return false
        c2 = size(file["energies"],2) ==12
        return c2
    end
end
function getRes(folder,L)
    files = let
        filesunsrt = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
        filter!(contains("Random_shuffle"),filesunsrt)
        filter!(contains("L=$L"),filesunsrt)
        # validfiles = is_valid_file.(filesunsrt)
        # if !all(validfiles)
        #     println("invalid files:")
        #     println(filesunsrt[.!validfiles])
        # end
        filter!(is_valid_file,filesunsrt)
        mus = [h5read(file,"mu") for file in filesunsrt]
        filesunsrt[sortperm(mus)]
    end
    energies = [h5read(file,"energies") for file in files]
    energies = [e[1:minimum(x->size(x,1),energies),:] for e in energies]
    energies = stack(energies)
    mus = [h5read(file,"mu") for file in files]
    Sqs = stack([h5read(file,"SqsGFMC") for file in files])
    taus = [h5read(file,"tau") for file in files]
    p_Sq = stack([h5read(file,"p_Sq") for file in files])
    return (;energies,mus,Sqs,taus,files,p_Sq)
end
KPoints = Dict([
    "Γ" => SVector(0,0),
    "X" => SVector(pi,0),
    "M" => SVector(pi,pi),
    "X'" => SVector(0,pi)
    ])
##
outfilesFolder = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/"

res = Dict(L=>getRes(outfilesFolder,L) for L in (20,24,28))

##
function filterRes(res,key)

    inds = findall(contains(key),res.files)

    return @views (;energies = res.energies[:,:,inds],mus = res.mus[inds],Sqs = res.Sqs[:,:,:,inds],taus = res.taus[inds],files = res.files[inds])
end

with_theme(theme_SimpleTicks()) do
    ind = 1
    L = 24
    Nsites = length(res[L].Sqs[1:end-1,1:end-1,1,ind,begin])
    enmean = mean(res[L].energies,dims=2)[1:50,1,ind] ./ Nsites

    enstd = std(res[L].energies,dims=2)[1:50,1,ind] ./ Nsites
    # enmean = mean(entest,dims=2)[:,1]
    # enstd = std(entest,dims=2)[:,1]
    tau = res[L].taus[ind]
    errlines((eachindex(enmean).-1) .*tau,enmean,enstd,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # lines((eachindex(enmean).-1) .*tau,enmean,axis = (;ylabel = L"E/N_\text{sites}",xlabel = L"\tau"))
    # band!((eachindex(enmean).-1) .*tau,enmean - enstd , enmean + enstd,color = (:black,0.2))
    # current_figure()

end


##

with_theme(theme_SimpleTicks()) do
    L = 20
    muIndex = findfirst(>=(0.5),res[L].mus)
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
    muIndex = findfirst(>=(0.8),res[L].mus)
    SqsGFMC = res[L].Sqs[:,:,100,:,muIndex]./ 4
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
    return fig
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

    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    # Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end
