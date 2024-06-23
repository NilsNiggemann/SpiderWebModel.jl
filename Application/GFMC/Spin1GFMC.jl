#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=S1GFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
L = 40
nBra = 8
NSteps = 120_000
equilibration_steps = 50_000
NWalkers = 48*100
w_avg_estimate = 9*L #estimated average weight for each iteration, to reduce floating point errors
λ = 1
outfile = "/p/scratch/pmfrg/niggemann1/Spiderweb/Data3/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_Lam=$(λ)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##

parentState = SW.stencilConfig(zeros(L,L),1;
# boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)

gfuncparams = h5read("/p/scratch/pmfrg/niggemann1/Spiderweb/DataStochRec_open/StochRec_L=40_nBra=8_NW=4800_Lam=1.h5","2/params")
ψG = SW.LocalPlaquetteGuidingFunction(gfuncparams)
##
@info "starting run"  
@time results = SW.startManyWalkerGFMC(parentState,NWalkers,NSteps,nBra,ψG,λ;outfile,equilibration_steps,w_avg_estimate,pre_equilibration_steps=60*equilibration_steps)
