#!/bin/bash
#=
#!/bin/bash

#SBATCH --job-name=S1GFMC                 # replace name
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=3-00:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMC/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)


i_arg = parse(Int, ARGS[1])


L = 40

infile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic/StochRec_L=40_nBra=12_NW=6400.h5"
##
NSteps = 15_000
equilibration_steps = 2_000
NWalkers = 128*80
λ = 0
nBra = 14
scatter_fraction = 0.7
outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_2/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##
optim_params = h5read(infile,"params_steps")[:,:,23]
@assert !iszero(optim_params)
ψG = SW.FullVariationalGuidingFunction(optim_params)
##

parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
w_avg_estimate = 0.266*prod(size(parentState))
DT = SW.DiscreteTimeMethod(0,nBra,w_avg_estimate)

##
@info "starting run" L λ nBra NSteps NWalkers outfile

flush(stdout)
@time results = SW.startManyWalkerGFMC(parentState,DT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=500_000,scatter_fraction,outfile)