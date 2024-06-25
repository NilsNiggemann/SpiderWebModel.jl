#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=GFMCS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=3-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/SFW/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/SWF_Spin1GFMC_periodic.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

infile = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small/") for file in files][i_arg]

##

##
nBra = h5read(infile,"nBranch")
SaveConfigs = SW.readMMapArray(infile,"SaveConfigs");
L,_,NWalkers,NSteps = size(SaveConfigs)
λ = h5read(infile,"Λ")

outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_ForwardWalking/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##
ψG = SW.FullVariationalGuidingFunction(h5read(infile,"FullVariationalGuidingFunction/params"))
##

parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
w_avg_estimate = SW.h5read(infile,"w_avg_estimate")
DT = SW.DiscreteTimeMethod(λ,nBra,w_avg_estimate)

##
@info "starting run" L λ nBra NSteps NWalkers outfile

@time results = SW.startManyWalkerGFMC(parentState,DT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=100_000,scatter_fraction,outfile)