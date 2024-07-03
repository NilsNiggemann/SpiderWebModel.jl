#!/bin/bash
#=
#!/bin/bash
#SBATCH --job-name=FWS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=2-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/SFW/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=END
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/B2CorrForwardWalking_incl.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

L = 20

infiles = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic/") for file in files]

infile = only(filter(contains("L=$L"),infiles))
##
nBra = 8
NSteps = 15_000
equilibration_steps = 15_000
NWalkers = 128*4
λ = 0
scatter_fraction = 0.8
pre_equilibration_steps = 100_000
##
# testing params
# L = 8
# nBra = 1
# NSteps = 10
# equilibration_steps = 0
# NWalkers = 1
# λ = 0
# scatter_fraction = 0.8
# pre_equilibration_steps = 0
##
outfileGFMC = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small2/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
mkpath(dirname(outfileGFMC))
@assert !isfile(outfileGFMC) "file already exists!"

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
DT = SW.DiscreteTimeMethod(λ,nBra,w_avg_estimate)

##
@info "starting run" L λ nBra NSteps NWalkers outfileGFMC

@time res = SW.startManyWalkerGFMC(parentState,DT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile = outfileGFMC)

##
# ____________________StraightForwardWalking_____________________
@info "starting straight forward walking"
##
# binnumber,fileindex =  i_arg %NumBins +1  , i_arg ÷ NumBins +1
##
outfileSFW = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_ForwardWalking_2/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
mkpath(dirname(outfileSFW))
@assert !isfile(outfileSFW) "file already exists!"


#___________Spin-1_______________________
##

refPlaq = SW.getCentralPlaquette(parentState)
symReduc = SW.symmetryReducePlaquettes(parentState,refPlaq)
GFMCPlaqs = collect(SW.plaquetteIterator(parentState))[symReduc.uniqueInds]

mProj = 300 ÷ nBra
BOp = SW.PlaquetteFlipOperator(parentState)
BBOp = SW.BBOperator(parentState,refPlaq)
##
@info "starting run" L λ mProj nBra NSteps NWalkers outfileSFW

resB = SW.measure_operator(parentState,DT,res.SaveConfigs,mProj,BOp,ψG,[refPlaq];outfile = outfileSFW)
##
resBB = SW.measure_operator(parentState,DT,res.SaveConfigs,mProj,BBOp,ψG,GFMCPlaqs;outfile = outfileSFW)