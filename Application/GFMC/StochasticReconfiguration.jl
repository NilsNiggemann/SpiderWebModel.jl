#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=SRS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/SR/SR_1_%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0
~/.bashrc
export JULIA_DEPOT_PATH="$SCRATCH/.julia"
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/StochasticReconfiguration.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
L = [40,30,28,24,20][i_arg]
nBra = 8
equilibration_steps = 1_000
Nwalkers = 128*7
λ = 0
outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic/L=($L)/StochRec_L=$(L)_nBra=$(nBra)_NW=$(Nwalkers).h5"

mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##
parentState = SW.stencilConfig(zeros(L,L),1,
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
ψG = SW.fullVariationalFunction(parentState,0.15)
DT = SW.DiscreteTimeMethod(0.,nBra,prod(size(parentState))*0.262)

##
@info "starting run"  
NSteps = 800
NBins = 30

stochReconfRes = SW.stochastic_reconfiguration(parentState,DT,i->min(NSteps + 200*i,5000),ψG,NBins,i -> max(0.8 - 0.1log(i),0.2),SW.IterativeSRSolver();Nwalkers,rel_tolerance=1e-8,equilibration_steps,pre_equilibration_steps=50_000,scatter_fraction=0.5,outfile)
