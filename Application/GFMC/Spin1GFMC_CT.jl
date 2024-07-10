#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8712861
#SBATCH --job-name=CTS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=3-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCT/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC_CT.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

Ls = [12,14,16,20]
L = Ls[(i_arg-1)%length(Ls)+1]
infiles = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic/") for file in files]
infile = only(filter(contains("L=$L"),infiles))
##
τ = 0.2
NSteps = 50_000
equilibration_steps = 5_000
NWalkers = 128*200
scatter_fraction = 0.8
outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"

mkpath(dirname(outfile))

@assert !isfile(outfile) "file already exists!"


##

using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##

optim_params = h5read(infile,"params_steps")[:,:,end]
@assert !iszero(optim_params)
ψG = SW.FullVariationalGuidingFunction(optim_params)
##
parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
w_avg_estimate = 0.266*prod(size(parentState))
CT = SW.ContinuousTimeMethod(τ,1,w_avg_estimate)


##
@info "starting run" L τ nBra NSteps NWalkers outfile

@time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=100_000,scatter_fraction,outfile)