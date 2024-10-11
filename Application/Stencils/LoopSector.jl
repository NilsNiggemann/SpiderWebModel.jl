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
    S_cond .= SW.periodicStateLoops(size(S,1))
    return S_cond
end
function get_S_diag(S)
    S_diag = copy(S)
    S_diag .= SW.periodicStateDiag(size(S,1))
    return S_diag
end
function get_S_stair(S)
    S_staircase = copy(S)
    S_staircase .= 4*SW.getStairCase(size(S,1))
    return S_staircase
end
function get_S_plainWeave(S)
    S_plainWeave = copy(S)
    S_plainWeave .= SW.periodicPlainWeave(size(S,1))
    return S_plainWeave
end


##
L = 20
println("Starting")
S = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
##
S_string = copy(S)
# S .= SW.h5read("temp.h5","conf")
S_string[end÷2,1:2:end] .= 2


S_two_strings = copy(S)
S_two_strings[end÷2,2:2:end] .= 2
S_two_strings[1:2:end,end÷2+1] .= -2


S_four_strings = copy(S)
S_four_strings[end÷4,1:2:end] .= 2
S_four_strings[2:2:end,end÷4] .= -2
S_four_strings[end÷4*3+1,1:2:end] .= 2
S_four_strings[1:2:end,end÷4*3] .= -2

# ψG = SW.fullVariationalFunction(S,0.15*(1-μ))

##
S_string_condensate = get_S_condensate(S)
S_diag_condensate = get_S_diag(S)
S_plainWeave = get_S_plainWeave(S)
S_stair = get_S_stair(S)
S_LoopsDense = copy(S) .= SW.periodicStateDenseLoops(size(S,1))

##
@assert SW.fulFillsConstraint(S_two_strings)
@assert SW.fulFillsConstraint(S_string_condensate)
@assert SW.fulFillsConstraint(S_diag_condensate)
@assert SW.fulFillsConstraint(S_plainWeave)
@assert SW.fulFillsConstraint(S_LoopsDense)
##
μ = 0.0
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))

##
function makeRuns(S_vec::AbstractVector,muRange,folder;Nwalkers = 30,NSteps = 8000,NwalkersOpt = 1,NStepsOpt = NSteps,OptIndep = Threads.nthreads())
    initializer0 = SW.WeightedConfigsInitializers(S_vec,[1 for _ in S_vec])
    S = first(S_vec)
    for μ in muRange
        ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
        
        files_folder = joinpath(folder,"μ=$(μ)")
        
        if !isdir(files_folder)
            mkpath(files_folder)
            println("Making runs for μ = $μ")
            flush(stdout)
            CTFindOpt = SW.ContinuousTimeMethod(1. + 2μ,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))
            
            @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NwalkersOpt,NStepsOpt,ψG;equilibration_steps=1,pre_equilibration_steps=40_000,scatter_fraction=0.5,initializer = initializer0) for _ in 1:OptIndep])
            initializer = SW.WeightedConfigsInitializers(OptimStart)
            flush(stdout)
            CT2 = SW.ContinuousTimeMethod(0.08,1,0.266*length(S)*(1-μ),SW.Hxx_RK(μ))
            
            @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,Nwalkers,NSteps,ψG;equilibration_steps=0,initializer=initializer,outfile = joinpath(files_folder,"$i.h5")) for i in 1:6])
            
            # plotEnergies(results,CT2;normalize=true,dense=true,τ = 10,color = :red,legend = false,axis = (;title = L"μ = %$μ"))
            # display(current_figure())
            # GC.gc()
            flush(stdout)
        end
    end
end
makeRuns(S::SW.StencilSpinConfig,args...;kwargs...) = makeRuns([S],args...;kwargs...)

# muRange = [0,0.02,0.025,0.03,0.035,0.1,0.15,0.2,0.3]
# muRange = [0.0,0.1,0.02,0.035,0.05,0.06,0.07,0.08,0.9,0.1,0.15]
# muRange = collect(0.0:0.02:0.12)
# muRange = collect(0.0:0.05:1)
muRange = collect(-0.1:0.05:1)
# append!(muRange,collect(0.05:0.01:0.09))
sort!(muRange,rev = true)
# folder = "temp/LoopSector/L=$(size(S,1))/noString"
# S_hash = get_S_hash(S)
# makeRuns([S,S_hash],muRange,folder;Nwalkers = 200,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 300,OptIndep = 12)
##
# folder = "temp/LoopSectorHybridInit2/L=$(size(S,1))/string"
# makeRuns(S_string,muRange,folder;Nwalkers = 80,NSteps = 2000,NwalkersOpt = 120,NStepsOpt = 100,OptIndep = 4)
# ##
# folder = "temp/LoopSector/L=$(size(S,1))/two_strings"
# makeRuns(S_two_strings,muRange,folder;Nwalkers = 180,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 4)
##
# folder = "temp/LoopSector/L=$(size(S,1))/four_strings"
# makeRuns(S_four_strings,muRange,folder;Nwalkers = 400,NSteps = 4000,NwalkersOpt = 600,NStepsOpt = 1000,OptIndep = 24)
##
# folder = "temp/LoopSectorNoInit/L=$(size(S,1))/string_condensate"
# makeRuns(S_string_condensate,muRange,folder;Nwalkers = 120,NSteps = 4000,NwalkersOpt = 160,NStepsOpt = 100,OptIndep = 4)
##
# folder = "temp/LoopSector/L=$(size(S,1))/diag_condensate"
# makeRuns(S_diag_condensate,muRange,folder;Nwalkers = 180,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 4)

##
# folder = "temp/LoopSector/L=$(size(S,1))/stair"
# makeRuns(S_stair,muRange,folder;Nwalkers = 500,NSteps = 10000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 6)
##
# folder = "temp/LoopSector/L=$(size(S,1))/plainWeave"
# makeRuns(S_plainWeave,muRange,folder;Nwalkers = 150,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 4)
##
folder = "temp/LoopSector/L=$(size(S,1))/LoopsDense"
makeRuns(S_LoopsDense,muRange,folder;Nwalkers = 150,NSteps = 4000,NwalkersOpt = 360,NStepsOpt = 100,OptIndep = 4)

##
if "TERM_PROGRAM" ∉ keys(ENV)
    exit()
end
##
# import Interpolations as ITP

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

res_noString_20 = getMuSweep("temp/LoopSector/L=20/noString/",100)
res_two_strings_20 = getMuSweep("temp/LoopSector/L=20/two_strings",100)
res_string_condensate_20 = getMuSweep("temp/LoopSector/L=20/string_condensate",100)
res_diag_condensate_20 = getMuSweep("temp/LoopSector/L=20/diag_condensate",100)
res_stair_20 = getMuSweep("temp/LoopSector/L=20/stair",100)
res_plainWeave_20 = getMuSweep("temp/LoopSector/L=20/plainWeave",100)
# res_plainWeave_40 = getMuSweep("temp/LoopSector/L=40/plainWeave",100)
res_LoopsDense_20 = getMuSweep("temp/LoopSector/L=20/LoopsDense",100)
##
function EnergyPlot(resultsVec,labels,L;
    colors = Makie.Colors.distinguishable_colors(length(resultsVec),parse.(Makie.Colors.RGB,[:black,:red,:blue])),
    markersize = 7,
    linestyle  = :solid,
    marker = '●'
    )
    QCP = 0.25
    with_theme(theme_SimpleTicks()) do 

        # res_four_strings = getMuSweep("temp/LoopSector3/L=14/four_strings")
        # res_string_condensate = getMuSweep("temp/LoopSector3/L=14/string_condensate")

        fig = Figure(size = 150 .* (4,5))
        MU_SCALE = 0.266
        comparison_func(μ) = MU_SCALE * (μ-1)
        
        xticks = collect(0:0.5:1)

        push!(xticks,QCP)
        unique!(sort!(xticks))
        xticklabels = [x == QCP ? L"μ_c = %$QCP" : L"%$x" for x in xticks]

        xminorticks = collect(0:0.1:1)
        ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"E_0/L^2",xlabelvisible=false,xticklabelsvisible=false,xminorticksvisible = true,xminorticks = xminorticks)

        axDiff = Axis(fig[2,1],xlabel = L"μ",ylabel = L"E_0/L^2 - %$MU_SCALE(μ-1)",xminorticks = xminorticks,xminorticksvisible = true,xticks = (xticks,xticklabels),
        # xlabelvisible=false,xticklabelsvisible=false
        )

        axderiv= Axis(
            fig[1,1],xlabel = L"μ",ylabel = L"dE/dμ",xminorticks = xminorticks,xminorticksvisible = true,xticks = (xticks,xticklabels),
            yaxisposition=:right,yticklabelcolor=:red,
            yticks = SimpleTicks(),
            xlabelvisible = false,
            ygridvisible=false,
            xgridvisible=false,
            xticklabelsvisible= false,
            xticksvisible= false,
            # xminorticksvisible=false
            ylabelvisible=true,
            ylabelcolor = :red
        )
        linkxaxes!(ax,axDiff,axderiv)
        Nsites = 20^2
        
        function getdEDmu(mu,E)
            # E_interp_linear = ITP.interpolate((mu,), E, ITP.Gridded(ITP.Linear()))
            # x -> ITP.gradient(E_interp_linear,x)[1]
            diff(E) ./ diff(mu)
        end

        function getDeltaDE(Estd)
            ΔDE = [Estd[i] + Estd[i+1] for i in eachindex(Estd)[1:end-1]]
        end

        # for (L,res_s,linestyle,marker,markersize) in zip(Ls,[res_s_20],linestyles,markers,markersizes)
        Nsites = L^2
        # getE(res) = minimum(res.E,dims=2)[:] ./ Nsites
        getEScale(E) = @views E[end,:] ./ Nsites
        getEDiff(res) = getEScale(res.E) .- comparison_func.(res.mus)

        for (res,color,label) in zip(resultsVec,colors,labels)

            E_mu = getEScale(res.E)
            E_std = getEScale(res.Estd)
            E_diff = getEDiff(res)
            μ = res.mus

            scatterlines!(ax,μ,E_mu;color ,label,markersize,linestyle,marker = marker)
            errorbars!(ax,μ,E_mu, E_std,color = color,whiskerwidth = 6,linewidth=0.5)

            scatterlines!(axDiff,μ,E_diff;color ,label,markersize,linestyle,marker = marker)
            errorbars!(axDiff,μ,E_diff, E_std,color = color,whiskerwidth = 6,linewidth=0.5)
            
            mu_fine = LinRange(extrema(μ)...,100)
            # dE = getdEDmu(res.mus,E_mu)
            # dE_interp = dE.(mu_fine)
            # μ = μ[1:2:end]
            # E_mu = E_mu[1:2:end]
            # E_std = E_std[1:2:end] 

            dE = getdEDmu(μ,E_mu)
            deltaDE = getDeltaDE(E_std./Nsites)

            scatterlines!(axderiv,μ[1:end-1],dE;color ,label,marker = marker,markersize = markersize,linestyle)
            errorbars!(axderiv,μ[1:end-1],dE,deltaDE,color = color,alpha = 0.2,whiskerwidth = 6,linewidth=0.5)
            # band!(axderiv,res.mus[begin:end-1],dE,deltaDE,color = color,alpha = 0.2)
            # lines!(ax,[0,1],[res.E[end,1] /Nsites,0],color = color,linestyle = :dash)
        end
        scatterlines!(ax,[NaN],[NaN];marker,linestyle,label = L"L = %$L",color = :grey)
        axislegend(ax,position = :rc,merge = true,unique = true,nbanks=1)
        vlines!(ax,[QCP],color = :grey,linestyle = :dash)
        vlines!(axDiff,[QCP],color = :grey,linestyle = :dash)
        rowsize!(fig.layout,1,Relative(0.6))
        fig
        
    end
        # xlims!(axDiff,-0.005,0.15)
    # ylims!(ax,-0.275,-0.22)

    
end
##
let 
    # colors = [:black,:red,:blue,:green,:orange,:purple]
    labels = [
        L"condensate$$",
        L"ℓ=0",
    # L"ℓ = 1",
        # L"ℓ = 2",
        L"diag$$",
        L"stair$$",
        L"plainWeave$$",
        L"LoopsDense$$",
    # "ℓ = 4",
    ]
    EnergyPlot([
        res_string_condensate_20,
        res_noString_20,
        # res_two_strings_20,
        res_diag_condensate_20,
        res_stair_20,
        res_plainWeave_20,
        res_LoopsDense_20
    ],labels,20)
end
##
EnergyPlot([res_plainWeave_40],[L"1"],40)


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
μ = 0.2
# res1 = getRes("temp/LoopSector/L=20/noString/μ=$μ")
# res1 = getRes("temp/LoopSector/L=20/string_condensate/μ=$μ")
# res1 = getRes("temp/LoopSector/L=20/diag_condensate/μ=$μ",3600)
# res1 = getRes("temp/LoopSector/L=20/stair/μ=$μ")
# res1 = getRes("temp/LoopSector/L=20/string_condensate/μ=$μ")
res1 = getRes("temp/LoopSector/L=20/LoopsDense/μ=$μ")
# res1 = getRes("temp/LoopSector/L=40/plainWeave/μ=$μ")
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
Snew = SW.stencilConfig(zeros(24,24),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
# Snew[end÷4,2:2:end] .= -2
# Snew[2:2:end,end÷4] .= -2
# Snew[end÷4*3,2:2:end] .= 2
# Snew[1:2:end,end÷4*3] .= 2
SW.flipSpinsAlongDiagonal!(Snew,23,-1;add=2)
Snew[24,23] = 2
Snew[23,24] = 2
# @assert SW.fulFillsConstraint(Snew)
SW.plotFractons(Snew)
##
μ = 0.1
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
CT = SW.ContinuousTimeMethod(0.08,1,(1-μ)* 0.9*0.266*length(Snew),SW.Hxx_RK(μ))
@time results = fetch.([SW.startManyWalkerGFMC(Snew,CT,120,1000,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5) for _ in 1:2])

##
plotEnergies(results,CT;normalize=true,dense=true,τ = 10)
##
SqsGFMC = mean(SW.getSqsGFMC(results,20,nBra = 1))

with_theme(theme_PiTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y")
    Sq = SW.getSqCont(SqsGFMC)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqsGFMC,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    fig
end
##
S = SW.stencilConfig(zeros(18,18),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)


S_hash = get_S_hash(S)
# @assert SW.fulFillsConstraint(S_hash)
μ = 0.2

SW.plotFractons(S_hash)
##

# CTFindOpt = SW.ContinuousTimeMethod(1.,1,(1-μ)* 0.9*0.266*length(S),SW.Hxx_RK(μ))
            
# @time OptimStart = fetch.([SW.startManyWalkerGFMC(S,CTFindOpt,200,100,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5) for _ in 1:12])
initializer = SW.WeightedConfigsInitializers([S,S_hash],[1,1])

##
CT2 = SW.ContinuousTimeMethod(0.08,1,0.266*length(S)*(μ-1),SW.Hxx_RK(μ))
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
res1 = fetch.([SW.startManyWalkerGFMC(S_hash,CT2,50,2000,ψG;equilibration_steps=0,initializer) for _ in 1:3])
res2 = fetch.([SW.startManyWalkerGFMC(S,CT2,50,2000,ψG;equilibration_steps=0,pre_equilibration_steps = 10_000,scatter_fraction = 0.0) for _ in 1:3])
##

plotEnergies(res1,CT2;normalize=true,dense=true,τ = 10,color = :green,nThermal = 100)
plotEnergies!(res2,CT2;normalize=true,dense=true,τ = 10,color = :red)
# plotEnergies(vcat(res1,res2),CT2;normalize=true,dense=true,τ = 10,color = :black,nThermal = 50)
current_figure()
##
equilib_plots(res1;scatter_fraction=0.8)
# equilib_plots(hcat(res1,res2);scatter_fraction=0.8)
##
SqsGFMC = mean(SW.getSqsGFMC(vcat(res1,res2),80,nBra = 1))

##
with_theme(theme_PiTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y")
    Sq = SW.getSqCont(SqsGFMC)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqsGFMC,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    fig
end
##
function round_matrix_elements(array, tol = 0.1)
    scale = maximum(abs,array)

    
    for i in eachindex(array)
        element = array[i]
        if abs(element) < tol*scale
            array[i] = 0
        else
            array[i] = sign(element)
        end
    end
    
    return array
end