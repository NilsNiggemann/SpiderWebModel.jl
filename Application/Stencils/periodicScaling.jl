cd(@__DIR__)
using Pkg
Pkg.activate("../.")
using Statistics, HDF5
import SpiderWebModel as SW
##
files = readdir("/storage/niggeni/Spiderweb/Data2/L=16/",join=true)
filesL = Dict(L => filter(contains("L=$L"),files) for L in [16,])
##
AllResults = Dict(L=> vcat(SW.readResults.(files,50_000)...) for (L,files) in filesL);
##
function getEns(results)
    en = stack([SW.getEnergies(res.TotalWeights,res.energies,1,550÷res.nBra) for res in results])
end
ens = Dict(L => getEns(results) for (L,results) in AllResults)

outfile = "../../Data2/Spin1GFMC_Eval_scaling.h5"
mkpath(dirname(outfile))
h5open(outfile,"w") do file
    for (L,en) in ens
        file["L=$L/energies"] = en
    end
end
##

let 
    for (L,res) in AllResults
        # for projectionSteps in (100,)
        for projectionSteps in (2000,1000,500)
            println((;L,projectionSteps))
            @time SqsGFMC = stack(SW.getSqsGFMC(res,projectionSteps),dims=3)
            h5open(outfile,"cw") do file
                file["SqsGFMC/L=$L/$projectionSteps"] = SqsGFMC
            end
        end
    end
end

##
using CairoMakie, MakieHelpers

outfile = "../../Data2/Spin1GFMC_Eval_scaling.h5"
let 
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1])
    ens = Dict(L=> h5read(outfile,"L=$L/energies") for L in [16,])
    normalize = true
    
    for (L,Result) in ens
        normalize || (L = 1)
        E = mean(Result,dims=2)[:] ./L.^2 
        lines!(ax,8 .*eachindex(E),E,label="L=$L")
        errorbars!(ax,8 .*eachindex(E),E, sqrt.(std(Result,dims=2)[:]) ./L.^2)
    end
    axislegend(ax)
    fig
end

##
with_theme(theme_PiTicks()) do 
    fig = Figure(size =(600,400))
    hmaxes = [Axis(fig[1, i],aspect = 1) for i in 1:3]
    linesaxes = [Axis(fig[2, i],aspect = 1,yticks = SimpleTicks()) for i in 1:3]
    for (j,L) in enumerate(keys(AllResults))
        SqsGFMC = h5read(outfile,"SqsGFMC/L=$L/500")
        err = sqrt.(var(SqsGFMC,dims=3))[:,:,1]
        SqsGFMC = mean(SqsGFMC,dims=3)[:,:,1]
        k = 2pi .* collect(0:L) ./L
        heatmap!(hmaxes[j],k,k,SqsGFMC)
        lines!(linesaxes[j],k,SqsGFMC[1,:],label="L=$L")
        errorbars!(linesaxes[j],k,SqsGFMC[1,:],2*err[1,:])
    end
    fig
end