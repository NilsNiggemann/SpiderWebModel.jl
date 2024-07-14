#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8745821
#SBATCH --job-name=CTRKS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=4-10:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC_CT_RK.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

μs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
μ = μs[(i_arg-1)%length(μs)+1]
L = 30
τ = 0.1
nBra = 1
NSteps = 2_000
equilibration_steps = 3_000
NWalkers = 128*16
scatter_fraction = 0.6

infiles = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic_RK_2/") for file in files]
filter!(contains("L=$L"),infiles)
filter!(contains("mu=$(μ)"),infiles)
infile = only(infiles)
##
outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg).h5"

mkpath(dirname(outfile))



##

using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##

optim_params = h5read(infile,"params_steps")[:,:,end]
@assert !iszero(optim_params)
# ψG = SW.FullVariationalGuidingFunction(optim_params)
ψG = SW.PlaquetteNumberGuidingFunction(0.1)
w_avg_estimate = -h5read(infile,"E0")[end]
parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))


##
if !isfile(outfile)

    @info "starting run" L τ nBra NSteps NWalkers outfile

    @time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=200_000,scatter_fraction,outfile)
end
results = SW.readResults(outfile,NSteps)[1]
##

@info "starting straight forward walking"
##
# binnumber,fileindex =  i_arg %NumBins +1  , i_arg ÷ NumBins +1
##
outfileSFW = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_SFW/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg).h5"
mkpath(dirname(outfileSFW))
@assert !isfile(outfileSFW) "file already exists!"


GFMCPlaqs = collect(SW.plaquetteIterator(parentState))

mProj = round(Int,15 ÷ τ)

BOp = SW.RandomPlaquetteFlipOperator(parentState)
##
@info "starting run" L τ  mProj nBra NSteps NWalkers outfileSFW

resB = SW.measure_operator(parentState,CT,results.SaveConfigs,mProj,BOp,ψG,[first(GFMCPlaqs)];outfile = outfileSFW)
##
Gnp = SW.precomputeNormalizedAccWeight(results.TotalWeights,1,mProj)

BVals = SW.get_observables_sfw(Gnp,resB[:,1,:]',SW.mean(results.TotalWeights)) ./length(collect(SW.plaquetteIterator(parentState)))

outfileSFW2 = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_SFW/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg)_result.h5"

h5write(outfileSFW2,"B",BVals)