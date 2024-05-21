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

function readResults(filename,binsize)
    energies_raw = h5read(filename,"energies")
    TotalWeights_raw = h5read(filename,"TotalWeights")
    reconfTable_raw = SW.readMMapArray(filename,"reconfigurationTable")
    SaveConfigs_raw = SW.readMMapArray(filename,"SaveConfigs")
    nBra = h5read(filename,"nBranch")
    getrange(i) = i*binsize+1:(i+1)*binsize
    @views function getRes(range)
        energies = energies_raw[range]
        TotalWeights = TotalWeights_raw[range]
        reconfTable = reconfTable_raw[:,range]
        SaveConfigs = SaveConfigs_raw[:,:,:,range]
        return (;energies,TotalWeights,SaveConfigs,reconfTable,nBra)
    end
    return [getRes(getrange(i)) for i in 0:length(energies_raw)÷binsize-1]
end
binsize=55_000
files = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/Data2/",join=true)[3:3]
AllResults = vcat(readResults.(files,binsize)...);
##
function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
end
en  = stack(getEns(AllResults))
outfile = "../Data/Spin1GFMC_Eval2.h5"
mkpath(dirname(outfile))
h5open(outfile,"w") do file
    file["energies"] = en
end
##

@views function getSq(res,p)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)
    Sq = similar(Conf, ComplexF64)
    
    Si = similar(Conf, ComplexF64)
    plan = SW.LatticeFFTs.FFTW.plan_fft(Conf)

    function SqFunc(Conf)
        Si .= Conf
        SW.mul!(Sq, plan, Si)
        Sq .= abs2.(Sq)
    end
    SaveConfs = res.SaveConfigs
    reconfTable = res.reconfTable
    res = SW.getObs(Gnp,SaveConfs,reconfTable,SqFunc,p÷2)
    newRes = similar(res,size(res).+1)
    newRes[begin:end-1,begin:end-1] .= res

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    newRes ./NSites
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getSqs(Results,p)
    Sqs = Vector{Matrix{Float64}}(undef,length(Results))
    Threads.@threads for i in eachindex(Results,Sqs)
        res = Results[i]
        nBra = res.nBra
        Sq = getSq(res,p÷nBra)
        Sqs[i] = Sq
    end
    return Sqs
end
for projectionSteps in (750,500,250)
# for projectionSteps in (20,40)
    SqsGFMC = stack(getSqs(AllResults,projectionSteps),dims=3)
    h5open(outfile,"cw") do file
        file["SqsGFMC/$projectionSteps"] = SqsGFMC
    end
end
##
