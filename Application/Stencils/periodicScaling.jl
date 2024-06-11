using Pkg
Pkg.activate("../.")
cd(@__DIR__)
# using CairoMakie, MakieHelpers
using Statistics, HDF5
import SpiderWebModel as SW
##
files = readdir("/storage/niggeni/Spiderweb/Data/",join=true)
filesL = Dict(L => filter(contains("L=$L"),files) for L in [16,18,20])
##
AllResults = Dict(L=> vcat(SW.readResults.(files,100_000)...) for (L,files) in filesL)
##
function getEns(results)
    en = stack([SW.getEnergies(res.TotalWeights,res.energies,1,250÷res.nBra) for res in results])
end
ens = Dict(L => getEns(results) for (L,results) in AllResults)

outfile = "../../Data/Spin1GFMC_Eval_scaling.h5"
mkpath(dirname(outfile))
h5open(outfile,"w") do file
    for (L,en) in ens
        file["L=$L/energies"] = en
    end
end
##
# let 
#     fig = Figure(size = (800, 600))
#     ax = Axis(fig[1, 1])
#     for (L,Result) in ens
#         E = mean(Result) ./L.^2
#         lines!(ax,8 .*eachindex(E),E,label="L=$L")
#         errorbars!(ax,8 .*eachindex(E),E, sqrt.(std(Result)) ./L.^2)
#     end
#     axislegend(ax)
#     fig
# end
##
let 
    for (L,res) in AllResults
        for projectionSteps in (50,100,200,250)
            println(L,projectionSteps)
            SqsGFMC = stack(SW.getSqsGFMC(res,projectionSteps),dims=3)
            h5open(outfile,"cw") do file
                file["SqsGFMC/L=$L/$projectionSteps"] = SqsGFMC
            end
        end
    end
end

