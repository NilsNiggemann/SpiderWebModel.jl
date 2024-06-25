#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=Spin1eval                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Spin1eval_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/Stencils/Spin1.jl ${SLURM_ARRAY_TASK_ID}
exit
=#
cd(@__DIR__)
using Pkg
Pkg.activate("../")

import SpiderWebModel as SW
using Statistics
using SpiderWebModel
using HDF5
##
# @views function readResults(filename,range)
#     energies = h5read(filename,"energies")[range]
#     TotalWeights = h5read(filename,"TotalWeights")[range]
#     reconfTable =  SW.readMMapArray(filename,"reconfigurationTable")[:,range]
#     SaveConfigs = SW.readMMapArray(filename,"SaveConfigs")[:,:,:,range]
#     nBra = h5read(filename,"nBranch")
#     return (;energies,TotalWeights,SaveConfigs,reconfTable,nBra)
# end


binsize=3_000

files = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small/L=28/") for file in files]

AllResults = vcat(SW.readResults.(files,binsize)...);
outfile = "../Data/Spin1GFMC_Eval_periodic.h5"
mkpath(dirname(outfile))
##

function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
end
en  = stack(getEns(AllResults))
##
h5open(outfile,"w") do file
    file["energies"] = en
    file["nBra"] = AllResults[1].nBra
    file["L"] = size(AllResults[1].SaveConfigs,1)
end
##

for projectionSteps in (1000,750,500)
# for projectionSteps in (20,40)
    SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
    h5open(outfile,"cw") do file
        file["SqsGFMC/$projectionSteps"] = SqsGFMC
    end
end
##
@views function getMagnetization(res,p)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = res.SaveConfigs[:,:,begin,begin]
    Si = zeros(size(Conf))

    function magFunc(Conf)
        Si .= Conf ./2
    end
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfTable
    res = SW.getObs(Gnp,SaveConfs,reconfTable,magFunc,p÷2)
end
function get_mags(Results,p)
    m = Vector{Matrix{Float64}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,m)
        res = Results[i]
        nBra = res.nBra
        Sq = getMagnetization(res,p÷nBra)
        m[i] = Sq
    end
    return m
end
##
for projectionSteps in (1000,750,500,250)
    m_GFMC = stack(get_mags(AllResults,projectionSteps),dims=3)
    h5open(outfile,"cw") do file
        file["magnetization/$projectionSteps"] = m_GFMC
    end
end