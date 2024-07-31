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

S = SW.stencilConfig(zeros(14,14),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

S_string = copy(S)
# S .= SW.h5read("temp.h5","conf")
S_string[end÷2,1:2:end] .= 2
# ψG = SW.fullVariationalFunction(S,0.15*(1-μ))
@assert SW.fulFillsConstraint(S_string)
# SW.plotFractons(S_string)
# SW.plotApplPlaquettes!(current_axis(),S_string)
# current_figure()
##

function makeRuns(S,muRange,folder)
    for μ in muRange
        ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
        
        files_folder = joinpath(folder,"μ=$(μ)")
        
        if !isdir(files_folder)
            mkpath(files_folder)
            CTFindOpt = SW.ContinuousTimeMethod(3.,1,(1-μ)* 0.2*0.266*length(S),SW.Hxx_RK(μ))
            
            @time OptimStart = SW.startManyWalkerGFMC(S,CTFindOpt,1200,1000,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5)
            
            initializer = SW.WeightedConfigsInitializers(OptimStart.SaveConfigs,OptimStart.TotalWeights)

            CT2 = SW.ContinuousTimeMethod(0.08,1,-mean(OptimStart.energies),SW.Hxx_RK(μ))
            
            @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,360,8000,ψG;equilibration_steps=0,initializer,outfile = joinpath(files_folder,"$i.h5")) for i in 1:6])
            
            # plotEnergies(results,CT2;normalize=true,dense=true,τ = 10,color = :red,legend = false,axis = (;title = L"μ = %$μ"))
            # display(current_figure())
            GC.gc()
        end
    end
end
muRange = [0,0.02,0.025,0.03,0.035,0.1,0.15,0.2,0.3]
append!(muRange,collect(0.05:0.01:0.09))
sort!(muRange)
folder = "temp/LoopSector3/L=$(size(S,1))/noString"
makeRuns(S,muRange,folder)
##
muRange = [0,0.02,0.04,0.05,0.15,0.2,0.3]
folder = "temp/LoopSector3/L=$(size(S,1))/string"
makeRuns(S_string,muRange,folder)
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
function getMuSweep(folder,N=100)
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
    folder = "temp/LoopSector3/L=14"
    res_noString = getMuSweep(folder*"/noString")
    res_string = getMuSweep(folder*"/string")
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"E_0/L^2")
    Nsites = length(S)
    scatterlines!(ax,res_noString.mus,res_noString.E[end,:] ./ Nsites,color = :red,label = L"no string$$")
    errorbars!(ax,res_noString.mus,res_noString.E[end,:] ./ Nsites,res_noString.Estd[end,:] ./ Nsites,color = :red,whiskerwidth = 10)
    scatterlines!(ax,res_string.mus,res_string.E[end,:] ./ Nsites,color = :blue,label = L"single string$$")
    errorbars!(ax,res_string.mus,res_string.E[end,:] ./ Nsites,res_string.Estd[end,:] ./ Nsites,color = :blue,whiskerwidth = 10)
    # band!(ax,res_noString.mus,res_noString.E[:,end] .- res_noString.Estd[:,end],res_noString.E[:,end] .+ res_noString.Estd[:,end],color = (:red,0.2))
    axislegend(ax,position = :rb)
    fig

    
end
