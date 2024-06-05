#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=B2S1GFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/B2Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/B2CorrForwardWalking.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)
##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5

i_arg = parse(Int, ARGS[1])

infiles = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/Data2/",join=true)

jobarray = [(file,range) for file in infiles for range in Iterators.partition(1:700_000,55_000)]
infile,Nrange = jobarray[i_arg]
@assert isfile(infile) "input file does not exist!"
##
nBra = h5read(infile,"nBranch")
mProj = 1000 ÷nBra
α = h5read(infile,"PlaquetteNumberGuidingFunction/alpha")
Λ = h5read(infile,"Λ")
w_avg_estimate = h5read(infile,"w_avg_estimate")
Lx,Ly,NWalkers,NSteps = h5open(infile,"r") do f
    f["SaveConfigs"] |> size
end

outfile = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataForwardWalking/Spin1GFMC_B2_L=$(Lx)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_alpha=$(α)_Lam=$(Λ)_$(i_arg).h5"
##
mkpath(dirname(outfile))
@assert !isfile(outfile) "output file already exists!"


#___________Spin-1_______________________
##
parentState = SW.stencilConfig(zeros(Lx,Ly),1)

GuidingWaveFunction = SW.PlaquetteNumberGuidingFunction(α)
##
@info "starting run" i_arg infile Nrange outfile Λ w_avg_estimate nBra mProj NWalkers NSteps length(jobarray)

# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,varFuncTest,1) for _ in 1:8])
# @time results = SW.startManyWalkerGFMC(parentState,NWalkers,NSteps÷nBra,nBra,GuidingWaveFunction,Λ)
SaveConfigs = @view SW.readMMapArray(infile,"SaveConfigs")[:,:,:,Nrange]
##
BOp = SW.PlaquetteFlipOperator(parentState)
resB = SW.measure_operator(parentState,SaveConfigs,mProj,nBra,BOp,GuidingWaveFunction,Λ,w_avg_estimate;outfile)
h5write(outfile,"NStart",first(Nrange))