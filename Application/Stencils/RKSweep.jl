using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
##
files = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK/") for file in files]
##
binsize=2_000

AllResults = stack(SW.readResults.(files[1:1],binsize));
    
##
function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
end
en  = stack(getEns(AllResults))
L = size(AllResults[1].SaveConfigs,1)
h5open(outfile,"cw") do file
    file["$L/energies"] = en
    file["$L/nBra"] = AllResults[1].nBra
    file["$L/L"] = L
end

    Threads.@threads for projectionSteps in (50,250,500,750,1000)
    # for projectionSteps in (20,40)
        println(L, " ",projectionSteps)
        SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
        h5open(outfile,"cw") do file
            file["$L/SqsGFMC/$projectionSteps"] = SqsGFMC
        end
    end
end
