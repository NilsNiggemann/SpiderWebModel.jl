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
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/B2CorrForwardWalking.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1]) - 1
NumBins = 2
binnumber,fileindex =  i_arg %NumBins +1  , i_arg ÷ NumBins +1

infile = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small/") for file in files][fileindex]

##

using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
##
nBra = h5read(infile,"nBranch")
L,_,NWalkers,NSteps = h5open(infile) do f
    size(f["SaveConfigs"])
end
res = SW.readResults(infile,NSteps÷ NumBins)[binnumber];
λ = h5read(infile,"Λ")
ψG = SW.FullVariationalGuidingFunction(h5read(infile,"FullVariationalGuidingFunction/params"))

outfile = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_ForwardWalking/L=$(L)/$i_arg/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"


#___________Spin-1_______________________
##

parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
w_avg_estimate = SW.h5read(infile,"w_avg_estimate")
DT = SW.DiscreteTimeMethod(λ,nBra,w_avg_estimate)

##

refPlaq = SW.getCentralPlaquette(parentState)
symReduc = SW.symmetryReducePlaquettes(parentState,refPlaq)
GFMCPlaqs = collect(SW.plaquetteIterator(parentState))[symReduc.uniqueInds]

mProj = 300 ÷ nBra
BOp = SW.PlaquetteFlipOperator(parentState)
BBOp = SW.BBOperator(parentState,refPlaq)
##
@info "starting run" L λ mProj nBra NSteps NWalkers outfile

resB = SW.measure_operator(parentState,DT,res.SaveConfigs,mProj,BOp,ψG,[refPlaq];outfile)
##
resBB = SW.measure_operator(parentState,DT,res.SaveConfigs,mProj,BBOp,ψG,GFMCPlaqs;outfile)