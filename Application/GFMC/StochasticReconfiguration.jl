#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=S1SR                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=mem192
#SBATCH --mem=170GB         # memory , more means less gc time
#SBATCH --time=0-12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Spin1StochRec_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/StochasticReconfiguration.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
L = [40,30,20][i_arg]
nBra = 8
equilibration_steps = 1_000
NWalkers = 48*100
λ = 1
outfile = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataStochRec_open/StochRec_L=$(L)_nBra=$(nBra)_NW=$(NWalkers)_Lam=$(λ).h5"

mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##
parentState = SW.stencilConfig(zeros(L,L),1)
ψG = SW.localPlaquetteGuidingFunction(parentState,0.13)
##
@info "starting run"  
NSteps = 2_000
NBins = 10

stochReconfRes = SW.stochastic_reconfiguration(parentState,NSteps,ψG,NBins,1e-3;Nwalkers = NWalkers,nbra = nBra,rel_tolerance=1e-1,equilibration_steps)

@info "step 1 completed"  

h5write(outfile,"1/params",stochReconfRes.params)
h5write(outfile,"1/E0_i",stochReconfRes.E0_i)
h5write(outfile,"1/ΔE_i",stochReconfRes.ΔE_i)
h5write(outfile,"1/params_steps",stochReconfRes.params_steps)
h5write(outfile,"1/NSteps",NSteps)
ψG = SW.LocalPlaquetteGuidingFunction(stochReconfRes.params)

##

NSteps = 5_000
NBins = 10

stochReconfRes = SW.stochastic_reconfiguration(parentState,NSteps,ψG,NBins,1e-2;Nwalkers = NWalkers,nbra = nBra,rel_tolerance=1e-1,equilibration_steps)

@info "step 2 completed"  

h5write(outfile,"2/params",stochReconfRes.params)
h5write(outfile,"2/E0_i",stochReconfRes.E0_i)
h5write(outfile,"2/ΔE_i",stochReconfRes.ΔE_i)
h5write(outfile,"2/params_steps",stochReconfRes.params_steps)
h5write(outfile,"2/NSteps",NSteps)

ψG = SW.LocalPlaquetteGuidingFunction(stochReconfRes.params)

##

NSteps = 8_000
equilibration_steps = 3_000
NBins = 30

stochReconfRes = SW.stochastic_reconfiguration(parentState,NSteps,ψG,NBins,4e-2;Nwalkers = NWalkers,nbra = nBra,rel_tolerance=1e-1,equilibration_steps)

@info "step 3 completed"  

h5write(outfile,"3/params",stochReconfRes.params)
h5write(outfile,"3/E0_i",stochReconfRes.E0_i)
h5write(outfile,"3/ΔE_i",stochReconfRes.ΔE_i)
h5write(outfile,"3/params_steps",stochReconfRes.params_steps)
h5write(outfile,"3/NSteps",NSteps)
