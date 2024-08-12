import Pkg
Pkg.activate(joinpath(@__DIR__,"../"))
cd(joinpath(@__DIR__,"../../"))
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
meanstd(x) = (mean(x),std(x))
function roundToEven(x)
    return Int(2*round(x/2))
end 

function get_S_hash(S)
    Base.require_one_based_indexing(S)
    S_hash = copy(S)
    Lx,Ly = size(S)
    
    S_hash[roundToEven(Lx÷4)+1,1:2:end] .= 2
    S_hash[roundToEven(Lx*3/4)+1,1:2:end] .= -2
    S_hash[2:2:end,roundToEven(Ly÷4)] .= 2
    S_hash[2:2:end,roundToEven(Ly*3/4)] .= -2
    return S_hash
end

function get_S_condensate(S)
    S_cond = copy(S)
    
    S_cond[begin:2:end,begin:4:end] .= 2

    S_cond[begin+1:4:end,begin+1:2:end] .= -2
    return S_cond
end

##
println("Starting")
S = SW.stencilConfig(zeros(22,22),1;
boundary = SW.Stencils.Remove(Int8(0)),padding = SW.Stencils.Halo(:in)
)
S_two_strings = copy(S)
S_two_strings[end÷2,2:2:end] .= 2
S_two_strings[1:2:end,end÷2+1] .= -2
S_string_condensate = get_S_condensate(S)

##
function makeRuns(S_vec::AbstractVector,muRange,folder;Nwalkers = 30,NSteps = 8000,NwalkersOpt = 1,NStepsOpt = NSteps,OptIndep = Threads.nthreads())
    initializer0 = SW.WeightedConfigsInitializers(S_vec,[1 for _ in S_vec])
    S = first(S_vec)
    for μ in muRange
        ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
        
        files_folder = joinpath(folder,"μ=$(μ)")
        
        if !isdir(files_folder)
            mkpath(files_folder)
            CTFindOpt = SW.ContinuousTimeMethod(1.,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))
            
            @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NwalkersOpt,NStepsOpt,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5,initializer = initializer0) for _ in 1:OptIndep])
            initializer = SW.WeightedConfigsInitializers(OptimStart)

            CT2 = SW.ContinuousTimeMethod(0.08,1,0.266*length(S)*(1-μ),SW.Hxx_RK(μ))
            
            @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,Nwalkers,NSteps,ψG;equilibration_steps=0,initializer=initializer,outfile = joinpath(files_folder,"$i.h5")) for i in 1:6])
            
            plotEnergies(results,CT2;normalize=true,dense=true,τ = 10,color = :red,legend = false,axis = (;title = L"μ = %$μ"))
            display(current_figure())
            GC.gc()
            flush(stdout)
        end
    end
end
makeRuns(S::SW.StencilSpinConfig,args...;kwargs...) = makeRuns([S],args...;kwargs...)

muRange = collect(0.0:0.05:1)
# append!(muRange,collect(0.05:0.01:0.09))
sort!(muRange)
folder = "temp/LoopSectorOpen/L=$(size(S,1))/noString"
S_hash = get_S_hash(S)
makeRuns([S,S_hash],muRange,folder;Nwalkers = 300,NSteps = 4000,NwalkersOpt = 400,NStepsOpt = 400,OptIndep = 16)
##
# folder = "temp/LoopSectorOpen/L=$(size(S,1))/two_strings"
# makeRuns(S_two_strings,muRange,folder;Nwalkers = 180,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 6)
##
# folder = "temp/LoopSectorOpen/L=$(size(S,1))/string_condensate"
# makeRuns(S_string_condensate,muRange,folder;Nwalkers = 120,NSteps = 4000,NwalkersOpt = 160,NStepsOpt = 100,OptIndep = 6)

##
# folder = "temp/LoopSectorNoInit/L=$(size(S,1))/string_condensate"
# makeRuns(S_string_condensate,muRange,folder;Nwalkers = 400,NSteps = 4000,NwalkersOpt = 600,NStepsOpt = 1000,OptIndep = 24)
if "TERM_PROGRAM" ∉ keys(ENV)
    exit()
end
##

function getAllFilesInFolder(folder)
    files = [joinpath(root,file) for (root,_,files) in walkdir(folder) for file in files]
    return files
end

function getRes(folder)
    files = getAllFilesInFolder(folder)
    results = [SW.readResults(file)[1] for file in files]
    return results
end

function getRes(folder,NMax)
    files = getAllFilesInFolder(folder)
    results = [SW.readResults(file,NMax)[1] for file in files]
    return results
end
function getMuSweep(folder,N=200,equilibration_steps=1)
    files = getAllFilesInFolder(folder)
    mus = sort!(unique([parse(Float64,match(r"μ=(-?\d+\.\d+)",file).captures[1]) for file in files]))
    E = zeros(N,length(mus))
    Estd = zeros(N,length(mus))
    for (i,mu) in enumerate(mus)
        filesmu = filter(contains("μ=$mu/"),files)
        results = [SW.readResults(file)[1] for file in filesmu]
        en = SW.getEnergies.(results,equilibration_steps,N)
        E[:,i] .= mean(en)
        Estd[:,i] .= std(en)
    end
    return (;mus,E,Estd)
end

# res_noString_18 = getMuSweep("temp/LoopSectorHybridInit/L=18/noString/")
# res_string_18 = getMuSweep("temp/LoopSectorHybridInit2/L=18/string")
# res_two_strings_18 = getMuSweep("temp/LoopSectorHybridInit2/L=18/two_strings")
# res_string_condensate_18 = getMuSweep("temp/LoopSectorNoInit/L=18/string_condensate")

res_noString_20 = getMuSweep("temp/LoopSectorOpen/L=20/noString/",50)
res_two_strings_20 = getMuSweep("temp/LoopSectorOpen/L=20/two_strings",50)
res_string_condensate_20 = getMuSweep("temp/LoopSectorOpen/L=20/string_condensate",50)
with_theme(theme_SimpleTicks()) do 

    # res_four_strings = getMuSweep("temp/LoopSector3/L=14/four_strings")
    # res_string_condensate = getMuSweep("temp/LoopSector3/L=14/string_condensate")

    fig = Figure(size = 150 .* (4,5))
    MU_SCALE = 0.266
    comparison_func(μ) = MU_SCALE * (μ-1)
    QCP = 0.25
    
    xticks = collect(0:0.5:1)

    push!(xticks,QCP)
    unique!(sort!(xticks))
    xticklabels = [x == QCP ? L"μ_c = %$QCP" : L"%$x" for x in xticks]

    xminorticks = collect(0:0.1:1)
    ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"E_0/L^2",xlabelvisible=false,xticklabelsvisible=false,xminorticksvisible = true,xminorticks = xminorticks)


    
    axDiff = Axis(fig[2,1],xlabel = L"μ",ylabel = L"E_0/L^2 - %$MU_SCALE(μ-1)",xminorticks = xminorticks,xminorticksvisible = true,xticks = (xticks,xticklabels))
    linkxaxes!(ax,axDiff)
    Nsites = 18^2
    
    colors = [:black,:red,:blue,:green,:orange]
    labels = [
        L"condensate$$",
        L"ℓ=0",
    # L"ℓ = 1",
        L"ℓ = 2",
    # "ℓ = 4",
    ]

    # res_s_18 = [
    #     res_string_condensate_18,
    #     res_noString_18,
    #     # res_string_18,
    #     res_two_strings_18,
    #     ]
        
    res_s_20 = [
            res_string_condensate_20,
            res_noString_20,
            res_two_strings_20,
    ]

    res_s_22 = empty(res_s_20)

    # getEDiff(res) = res.E[end,:] ./ Nsites
    linestyles = [:solid,:dash,:dot]
    markers = ('●','×' ,'+')
    Ls = [22,20]
    markersizes = [6,15]
    for (L,res_s,linestyle,marker,markersize) in zip(Ls,[res_s_22,res_s_20],linestyles,markers,markersizes)
        Nsites = L^2
        getEScale(res) = res.E[end,:] ./ Nsites
        getEDiff(res) = getEScale(res) .- comparison_func.(res.mus)

        for (res,color,label) in zip(res_s,colors,labels)
            scatterlines!(ax,res.mus,getEScale(res);color ,label,markersize,linestyle = linestyle,marker = marker)
            errorbars!(ax,res.mus,getEScale(res),res.Estd[end,:] ./ Nsites,color = color,whiskerwidth = 6,linewidth=0.5)

            scatterlines!(axDiff,res.mus,getEDiff(res);color ,label,markersize,linestyle = linestyle,marker = marker)
            errorbars!(axDiff,res.mus,getEDiff(res),res.Estd[end,:] ./ Nsites,color = color,whiskerwidth = 6,linewidth=0.5)
            
            # lines!(ax,[0,1],[res.E[end,1] /Nsites,0],color = color,linestyle = :dash)
        end
        scatterlines!(ax,[NaN],[NaN];marker,linestyle,label = L"L = %$L",color = :grey)

    end
    axislegend(ax,position = :lt,merge = true,unique = true,nbanks=1)
    # xlims!(axDiff,-0.005,0.15)
    # ylims!(ax,-0.275,-0.22)

    vlines!(ax,[QCP],color = :grey,linestyle = :dash)
    vlines!(axDiff,[QCP],color = :grey,linestyle = :dash)
    rowsize!(fig.layout,1,Relative(0.7))
    fig

    
end
##
using Optim
function SqFieldTheory(qx::Real, qy::Real, A, r)
    cos_qx = cos(qx)
    cos_qy = cos(qy)
    sin_qx = sin(qx)
    sin_qy = sin(qy)
    numerator = (2 * (cos_qx - cos_qy + 2 * sin_qx * sin_qy))^2
    denominator = sqrt(((cos_qx - cos_qy)^2 + 4 * sin_qx^2 * sin_qy^2) * (r + 4 * ((cos_qx - cos_qy)^2 + 4 * sin_qx^2 * sin_qy^2)))
    return A* numerator / (denominator+1e-30)
end

SqFieldTheory(q,v,w) = SqFieldTheory(q[1],q[2],v,w)

SqFieldTheory(q,coefs::AbstractVector) = SqFieldTheory(q[1],q[2],coefs[1],coefs[2])
function optimizeCoeffs(SqMat)
    q = trueMomenta(0., 2pi, size(SqMat, 1) - 1)[1:end-1]

    function loss(v, w)
        l = 0.0
        v = abs(v)
        w = abs(w)
        for (i, qx) in enumerate(q), (j, qy) in enumerate(q)
            l += abs2(SqMat[i, j] - SqFieldTheory(qx, qy, v, w))
        end
        return l
    end

    loss(v) = loss(v[1], v[2])

    x0 = [1., 1.]

    res = optimize(loss, x0)
    @info res
    coefs = abs.(Optim.minimizer(res))
    return coefs
end
function rasterCurve(curvePoints,grid,t)

    getPos(point) = findmin(x->SW.norm(SW.SVector(x.-point)),grid)[2]
    positions = getPos.(curvePoints)
    tnew = empty(t)
    posnew = empty(positions)
    for i in eachindex(t)
        p = positions[i]
        if p ∉ posnew
            push!(tnew, t[i])
            push!(posnew,p)
        end 
    end
    return tnew,posnew
end

using StaticArrays
function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,SW.norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end


##
# res1 = getRes("temp/LoopSector/L=20/noString/μ=0.0")
μ = 0.3
res1 = getRes("temp/LoopSectorOpen/L=20/noString/μ=$μ")
# res1 = getRes("temp/LoopSectorOpen/L=20/string_condensate/μ=$μ")
# res1 = getRes("temp/LoopSector/L=20/two_strings/μ=$μ")
##
SqsGFMC = SW.getSqsGFMC(res1,50)
SqMat = mean(SqsGFMC)
SqErr = std(SqsGFMC)
fittingCoefs = optimizeCoeffs(SqMat)
##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(size = 150 .* (4,5))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks(), yticks = PiTicks())
    axFT = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1,ylabelvisible = false,yticklabelsvisible = false,xticks = PiTicks(), yticks = PiTicks())
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
    qr = LinRange(0,.45pi,100)
    
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

    for (phi,color) in zip([0,pi/4,pi/2.8],colors)
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
    fig
end
##
Stest = SW.stencilConfig(zeros(20,20),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
Stest .= SW.periodicStateWeb(20)
SW.plotSpinConfig(Stest)