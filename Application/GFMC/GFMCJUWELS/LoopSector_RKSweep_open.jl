#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=OmuGFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-08:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK $EXECUTABLE_FILE ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

# muRange = [0,0.05,0.07,0.08,0.1,0.15,0.2]
muRange = collect(0:0.01:0.14)

μ = muRange[i_arg]
import Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using Statistics
using SpiderWebModel.HDF5
##
S = SW.stencilConfig(zeros(26,26),1)

# S_string = copy(S)
# # S .= SW.h5read("temp.h5","conf")
# S_string[end÷2,1:2:end] .= 2

# S_two_strings = copy(S)
# S_two_strings[end÷2,1:2:end] .= 2
# S_two_strings[2:2:end,end÷2+1] .= -2

# @assert SW.fulFillsConstraint(S_string)
# @assert SW.fulFillsConstraint(S_two_strings)
##

function makeRun(S,μ,folder;Nwalkers = 30,NSteps = 8000,NwalkersOpt = 1,NStepsOpt = NSteps,OptIndep = Threads.nthreads(),kwargs...)
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
    files_folder = joinpath(folder,"μ=$(μ)")
    
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
        
    files_folder = joinpath(folder,"μ=$(μ)")
    
    if !isdir(files_folder)
        mkpath(files_folder)
        CTFindOpt = SW.ContinuousTimeMethod(1.,1,(1-μ)* 0.22*length(S),SW.Hxx_RK(μ))
        
        @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NwalkersOpt,NStepsOpt,ψG;equilibration_steps=1,pre_equilibration_steps=30_000,scatter_fraction=0.5) for _ in 1:OptIndep])
        initializer = SW.WeightedConfigsInitializers(OptimStart)
        # initializer = initializer0

        CT2 = SW.ContinuousTimeMethod(0.1,1,(1-μ)* 0.22*length(S),SW.Hxx_RK(μ))
        
        @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT2,Nwalkers,NSteps,ψG;equilibration_steps=0,outfile = joinpath(files_folder,"$i.h5"),initializer,kwargs...) for i in 1:6])
        try
            println([mean(res.energies) for res in results])
        catch
            
        end
        flush(stdout)
        GC.gc()
    end
    GC.gc()
end

##
folder = ENV["MYSCRATCH"]*"/LoopSector_open/L=$(size(S,1))/noString"
makeRun([S,S_hash],μ,folder;Nwalkers = 48*70,NSteps = 5000,NwalkersOpt = 48*80,NStepsOpt = 400,OptIndep = 20,equilibration_steps=0)
GC.gc()