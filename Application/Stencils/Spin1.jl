#!/bin/bash
#=
#!/bin/bash


#SBATCH --job-name=Spin1eval                 # replace name
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=0-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/eval%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/Stencils/Spin1.jl ${SLURM_ARRAY_TASK_ID}
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

binsize=15000

groups = readdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_2/",join=true)

# outfile = "../Data/Spin1GFMC_Eval_periodic_L40.h5"
outfile = "/scratch/hpc-prf-pm2frg/niggeni/tmp/Data/Spin1GFMC_Eval_periodic_L40.h5"
mkpath(dirname(outfile))
function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
end
##
Threads.@threads for group in groups
    files = [joinpath(root,file) for (root,_,files) in walkdir(group) for file in files][1:1]
    @time AllResults = vcat(SW.readResults.(files,binsize)...);    


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
##
# @views function getMagnetization(res,p)
#     Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

#     Conf = res.SaveConfigs[:,:,begin,begin]
#     Si = zeros(size(Conf))

#     function magFunc(Conf)
#         Si .= Conf ./2
#     end
#     SaveConfs = res.SaveConfigs
#     reconfTable = res.reconfTable
#     res = SW.getObs(Gnp,SaveConfs,reconfTable,magFunc,p÷2)
# end
# function get_mags(Results,p)
#     m = Vector{Matrix{Float64}}(undef,length(Results))
#     Threads.@threads for i in eachindex(Results,m)
#         res = Results[i]
#         nBra = res.nBra
#         Sq = getMagnetization(res,p÷nBra)
#         m[i] = Sq
#     end
#     return m
# end
# ##
# for projectionSteps in (1000,750,500,250)
#     m_GFMC = stack(get_mags(AllResults,projectionSteps),dims=3)
#     h5open(outfile,"cw") do file
#         file["magnetization/$projectionSteps"] = m_GFMC
#     end
# end