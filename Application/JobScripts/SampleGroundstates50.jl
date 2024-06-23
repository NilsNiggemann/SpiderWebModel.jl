#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=GS50
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GS/GS50%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0; julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/JobScripts/SampleGroundstates50.jl $SLURM_ARRAY_TASK_ID
exit
=#
println(DEPOT_PATH)
cd(@__DIR__)
using Pkg
Pkg.activate(".")

import SpiderWebModel as SW
using SpiderWebModel.HDF5
##

function savefile(filename,stage,sols)
    h5open(filename,"cw") do f 
        f["Confs/$(stage)"] = sols
    end
end

function main()
    i_arg = parse(Int, ARGS[1])
    L = 50
    filename = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/Data/L_$L/Confs_$(L)_$i_arg.h5"
    mkpath(dirname(filename))
    TimeLimit = 100
    fixedFraction = 1/6
    NRuns = 10_000
    WorkLimit = 100
    for stage in 1:10
        newsols = stack(Matrix.(SW.constructGroundstates(L, NRuns, fixedFraction; TimeLimit,WorkLimit)),dims=3)
        savefile(filename,stage,newsols)
        @info "finished stage $stage" numSols = size(newsols,3)
    end
    return 0
end
main()

