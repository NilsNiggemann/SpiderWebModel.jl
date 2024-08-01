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
##

S = SW.stencilConfig(zeros(18,18),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

S_string = copy(S)
# S .= SW.h5read("temp.h5","conf")
S_string[end÷2,1:2:end] .= 2


S_two_strings = copy(S)
S_two_strings[end÷2,1:2:end] .= 2
S_two_strings[2:2:end,end÷2+1] .= -2

S_four_strings = copy(S)
S_four_strings[end÷4+1,1:2:end] .= 2
S_four_strings[2:2:end,end÷4] .= -2
S_four_strings[end÷4*3+1,1:2:end] .= 2
S_four_strings[2:2:end,end÷4*3] .= -2
##
S_string_condensate = copy(S)

# S .= SW.h5read("temp.h5","conf")
S_string_condensate[1:2:end,1:2:end] .= 2
S_string_condensate[2:2:end,2:2:end] .= -2
# ψG = SW.fullVariationalFunction(S,0.15*(1-μ))
##
@assert SW.fulFillsConstraint(S_string)
@assert SW.fulFillsConstraint(S_two_strings)
@assert SW.fulFillsConstraint(S_four_strings)
@assert SW.fulFillsConstraint(S_string_condensate)
##
SW.plotFractons(S_two_strings)
SW.plotApplPlaquettes!(current_axis(),S_two_strings)
current_figure()
##
μ = 0.1
ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))

CTFindOpt = SW.ContinuousTimeMethod(3.,1,(1-μ)* 0.9*0.266*length(S),SW.Hxx_RK(μ))
            
@time OptimStart = fetch.([SW.startManyWalkerGFMC(S,CTFindOpt,1,2000,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5) for _ in 1:6])
initializer = SW.WeightedConfigsInitializers(OptimStart,:energies)

##
function makeRuns(S,muRange,folder;Nwalkers = 30,NSteps = 8000,NwalkersOpt = 1,NStepsOpt = NSteps,OptIndep = Threads.nthreads())
    for μ in muRange
        ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
        
        files_folder = joinpath(folder,"μ=$(μ)")
        
        if !isdir(files_folder)
            mkpath(files_folder)
            CTFindOpt = SW.ContinuousTimeMethod(3.,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))
            
            @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NwalkersOpt,NStepsOpt,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5) for _ in 1:OptIndep])
            initializer = SW.WeightedConfigsInitializers(OptimStart)

            CT2 = SW.ContinuousTimeMethod(0.08,1,-mean(OptimStart[1].energies),SW.Hxx_RK(μ))
            
            @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,Nwalkers,NSteps,ψG;equilibration_steps=0,initializer,outfile = joinpath(files_folder,"$i.h5")) for i in 1:6])
            
            plotEnergies(results,CT2;normalize=true,dense=true,τ = 10,color = :red,legend = false,axis = (;title = L"μ = %$μ"))
            display(current_figure())
            GC.gc()
        end
    end
end
# muRange = [0,0.02,0.025,0.03,0.035,0.1,0.15,0.2,0.3]
muRange = [-0.1,0.0,0.1,0.02,0.035,0.05,0.06,0.07,0.08,0.9,0.1,0.15]
# append!(muRange,collect(0.05:0.01:0.09))
sort!(muRange)
folder = "temp/LoopSector/L=$(size(S,1))/noString"
makeRuns(S,muRange,folder;Nwalkers = 180,NSteps = 2000,NwalkersOpt = 600,NStepsOpt = 50,OptIndep = 12)
##
folder = "temp/LoopSector/L=$(size(S,1))/string"
makeRuns(S_string,muRange,folder)
##
folder = "temp/LoopSector/L=$(size(S,1))/two_strings"
makeRuns(S_two_strings,muRange,folder;Nwalkers = 180,NSteps = 3000,NwalkersOpt = 600,NStepsOpt = 50,OptIndep = 12)
##
folder = "temp/LoopSector/L=$(size(S,1))/four_strings"
makeRuns(S_four_strings,muRange,folder)

##
folder = "temp/LoopSector/L=$(size(S,1))/string_condensate"
makeRuns(S_string_condensate,muRange,folder)
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
function getMuSweep(folder,N=200)
    files = getAllFilesInFolder(folder)
    mus = sort!(unique([parse(Float64,match(r"μ=(\d+\.\d+)",file).captures[1]) for file in files]))
    E = zeros(N,length(mus))
    Estd = zeros(N,length(mus))
    for (i,mu) in enumerate(mus)
        filesmu = filter(contains("μ=$mu/"),files)
        results = [SW.readResults(file)[1] for file in filesmu]
        en = SW.getEnergies.(results,1,N)
        E[:,i] .= mean(en)
        Estd[:,i] .= std(en)
    end
    return (;mus,E,Estd)
end

with_theme(theme_SimpleTicks()) do 
    res_noString = getMuSweep("temp/LoopSector3/L=14/noString")
    res_string = getMuSweep("temp/LoopSector3/L=14/string")
    res_two_strings = getMuSweep("temp/LoopSector3/L=14/two_strings")
    # res_four_strings = getMuSweep("temp/LoopSector3/L=14/four_strings")
    res_string_condensate = getMuSweep("temp/LoopSector3/L=14/string_condensate")

    fig = Figure(size = 150 .* (4,3))
    ax = Axis(fig[1,1:4],xlabel = L"μ",ylabel = L"E_0/L^2")
    Nsites = length(S)
    
    colors = [:black,:blue,:red,:green,:orange]
    labels = [L"ℓ=0",
    L"ℓ = 1",
    L"ℓ = 2",
    # "ℓ = 4",
    L"ℓ=L^2"]

    res_s = [res_noString,
    res_string,
    res_two_strings,
    res_string_condensate]
    for (res,color,label) in zip(res_s,colors,labels)
        scatterlines!(ax,res.mus,res.E[end,:] ./ Nsites,color = color,label = label)
        errorbars!(ax,res.mus,res.E[end,:] ./ Nsites,res.Estd[end,:] ./ Nsites,color = color,whiskerwidth = 10)
        # lines!(ax,[0,1],[res.E[end,1] /Nsites,0],color = color,linestyle = :dash)
    end
    axislegend(ax,position = :rb)
    xlims!(ax,-0.005,0.15)
    ylims!(ax,-0.27,-0.2)

    axes = [Axis(fig[2,i];SW.getConfigAxis(S)...,title = labels[i],xticklabelsvisible = false,yticklabelsvisible = false) for (i,S) in enumerate([S,S_string,S_two_strings,S_string_condensate])]
    for (i,S) in enumerate([S,S_string,S_two_strings,S_string_condensate])
        SW.plotSpinConfig!(axes[i],S)
    end
    rowsize!(fig.layout,1,Relative(0.7))
    fig

    
end
##
res1 = getRes("temp/LoopSector3/L=14/noString/μ=0.15")
##
with_theme(theme_PiTicks()) do 
    SqsGFMC = mean(SW.getSqsGFMC(res1,200))
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
Snew = SW.stencilConfig(zeros(24,24),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
Snew[1,end] = -2
# Snew[end÷4,2:2:end] .= -2
# Snew[2:2:end,end÷4] .= -2
# Snew[end÷4*3,2:2:end] .= 2
# Snew[1:2:end,end÷4*3] .= 2
SW.flipSpinsAlongDiagonal!(Snew,23,-1;add=2)
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