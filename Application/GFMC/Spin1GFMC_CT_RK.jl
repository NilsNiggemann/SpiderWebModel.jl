#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8745821
#SBATCH --job-name=CTRKS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=3-10:00:00
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
τ = 0.2
nBra = 1
NSteps = 1_000
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

@assert !isfile(outfile) "file already exists!"


##

using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
#___________Spin-1_______________________
##

optim_params = h5read(infile,"params_steps")[:,:,1]
@assert !iszero(optim_params)
ψG = SW.FullVariationalGuidingFunction(optim_params)
# ψG = SW.PlaquetteNumberGuidingFunction(0.1)
w_avg_estimate = -h5read(infile,"E0")[end]
parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))


##
@info "starting run" L τ nBra NSteps NWalkers outfile

@time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=200_000,scatter_fraction,outfile)

##

@info "starting straight forward walking"
##
# binnumber,fileindex =  i_arg %NumBins +1  , i_arg ÷ NumBins +1
##
outfileSFW = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_SFW/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg).h5"
mkpath(dirname(outfileSFW))
@assert !isfile(outfileSFW) "file already exists!"


refPlaq = SW.getCentralPlaquette(parentState)
symReduc = SW.symmetryReducePlaquettes(parentState,refPlaq)
GFMCPlaqs = collect(SW.plaquetteIterator(parentState))[symReduc.uniqueInds]

mProj = round(Int,150 ÷ τ)

BOp = SW.PlaquetteFlipOperator(parentState)
##
@info "starting run" L λ mProj nBra NSteps NWalkers outfileSFW

resB = SW.measure_operator(parentState,DT,res.SaveConfigs,mProj,BOp,ψG,[refPlaq];outfile = outfileSFW)
##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj) for res in results]

BVals = stack([SW.get_observables_sfw(Gnp,res[:,1,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resB,results) ])

outfileSFW2 = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK_SFW/$i_arg/Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(i_arg)_result.h5"

h5write(outfileSFW2,"B",BVals)