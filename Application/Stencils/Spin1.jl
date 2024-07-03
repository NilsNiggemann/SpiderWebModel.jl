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

binsize=6_000

groups = readdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small",join=true)

outfile = "../Data/Spin1GFMC_Eval_periodic.h5"
mkpath(dirname(outfile))

Threads.@threads for group in groups
    files = [joinpath(root,file) for (root,_,files) in walkdir(group) for file in files]
    AllResults = vcat(SW.readResults.(files,binsize)...);
    

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

    Threads.@threads for projectionSteps in (50,100,200,500)
    # for projectionSteps in (20,40)
        println(L, " ",projectionSteps)
        SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
        h5open(outfile,"cw") do file
            file["$L/SqsGFMC/$projectionSteps"] = SqsGFMC
        end
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