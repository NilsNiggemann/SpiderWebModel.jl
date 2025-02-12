#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=S1GS40
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GS/S1GS40%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
# module --force purge
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/JobScripts/Spin1/SampleS1Groundstates40.jl $SLURM_ARRAY_TASK_ID
exit
=#
cd(@__DIR__)
using Pkg
Pkg.activate(".")

import SpiderWebModel as SW
using SpiderWebModel.HDF5
##

function savefile(filename,stage,sols,fixedFraction)
    h5open(filename,"cw") do f 
        f["Confs/$(stage)"] = sols
        f["fixedFraction/$(stage)"] = fixedFraction
    end
end

function main()
    i_arg = parse(Int, ARGS[1])
    L = 40
    filename = "$(ENV["MYSCRATCH"])/Spiderweb/Data/S1_L_$L/Confs_$(L)_$i_arg.h5"
    mkpath(dirname(filename))
    TimeLimit = 10000
    fixedFraction = 0.26
    NRuns = 8_000
    WorkLimit = 10000
    for stage in 21:40
        println("starting stage $stage")
        flush(stdout)
        fixedFraction += 0.005
        sols = SW.constructGroundstatesSpin1(L, NRuns, fixedFraction; STotZero=true,TimeLimit,WorkLimit)
        if !isempty(sols[1])
            newsols = stack(Matrix.(sols[1]),dims=3)
            savefile(filename,stage,newsols,fixedFraction)
        end
        @info "finished stage $stage" numSols = size(newsols,3)
        flush(stderr)
    end
    return 0
end
main()

