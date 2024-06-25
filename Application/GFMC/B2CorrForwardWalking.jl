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
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/B2CorrForwardWalking.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])

infile = [joinpath(root,file) for (root,_,files) in walkdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_small/") for file in files][i_arg]

##

using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
##
nBra = h5read(infile,"nBranch")
res = SW.readResults(infile,NSteps)[1];
L,_,NWalkers,NSteps = size(res.SaveConfigs)
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

##
using SpiderWebModel.Statistics
@info "run over. Starting evaluation step"
Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,mProj)

BBVals = [SW.get_observables_sfw(Gnp,resBB[:,j,:]',mean(res.TotalWeights)) for j in eachindex(GFMCPlaqs)]

BVals = SW.get_observables_sfw(Gnp,resBB[:,begin,:]',mean(res.TotalWeights))

##

using SpiderWebModel.StaticArrays

function getBBCorrelator(BBVals,BVals,symReduc;index = lastindex(BBVals[1]))

    BBCorrelatorRaw = getBBCorrelator(BBVals,BVals,index)

    BBCorrelator = similar(BBCorrelatorRaw, length(symReduc.indicesMapping))
    for (i,k) in enumerate(symReduc.indicesMapping)
        BBCorrelator[i] = BBCorrelatorRaw[k]
    end
    return BBCorrelator
end

function getBBCorrelator(BBVals,BVals,index::Int)
    nth(x) = x[index]
    BBEnd = nth.(BBVals)
    BEnd =  nth(BVals)
    BBCorrelatorRaw = BBEnd .- BEnd^2
    return BBCorrelatorRaw
end

BBCorrelators = stack([
    getBBCorrelator(BBVals,BVals,symReduc,index = i)
    for i in 1:mProj
    ]
)
##
outfileFinal = "../Data/BBCorr/BBCorr_L=$(L)_$(i_arg)"
mkpath(dirname(outfileFinal))
h5write(outfileFinal,"BBCorrelator",BBCorrelators)