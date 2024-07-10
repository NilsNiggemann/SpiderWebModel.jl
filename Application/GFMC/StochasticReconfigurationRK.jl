#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=SRS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=2-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/SR/SR_RK_%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0
~/.bashrc
export JULIA_DEPOT_PATH="$SCRATCH/.julia"
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/StochasticReconfigurationRK.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
μ = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9][i_arg]
L = 30
τ = 0.1
equilibration_steps = 1_000
Nwalkers = 128*5
λ = 0
outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic_RK_2/L=($L)/StochRec_L=$(L)_tau=$(τ)_NW=$(Nwalkers)_mu=$(μ).h5"

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
αstart = 0.14 * (1-μ)
ψG = SW.fullVariationalFunction(parentState,αstart)
CT = SW.ContinuousTimeMethod(τ,1,prod(size(parentState))*0.266*(1-μ),SW.Hxx_RK(μ))

##
@info "starting run"  
NSteps = 500
NBins = 200

stochReconfRes = SW.stochastic_reconfiguration(parentState,CT,i->min(NSteps + 50*i,3000),ψG,NBins,i -> min(0.5 + i,3.),SW.IterativeSRSolver();Nwalkers,rel_tolerance=1e-8,equilibration_steps,pre_equilibration_steps=50_000,scatter_fraction=0.5,outfile)
