#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=S1GFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)

i_arg = parse(Int, ARGS[1])
L = 40
nBra = 15
NSteps = 6_000_000
NWalkers = 48*10
α = 0.15
λ = 1

outfile = "/p/scratch/pmfrg/niggemann1/Spiderweb/Data2/Spin1GFMC_L=$(L)_nBra=$(nBra)_NSteps=$(NSteps)_NW=$(NWalkers)_alpha=$(α)_Lam=$(λ)_$(i_arg).h5"
mkpath(dirname(outfile))
@assert !isfile(outfile) "file already exists!"

##
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5, H5Zblosc
##
h5open(outfile,"w") do f
    f["test",blosc=9] = rand(Int8,20,20)
end
rm(outfile)
#___________Spin-1_______________________
##
parentState = SW.stencilConfig(zeros(L,L),1)

GuidingWaveFunction(x) = SW.varitationalFunc(α,x,0)
##
@info "starting run"  
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,varFuncTest,1) for _ in 1:8])
@time results = SW.startManyWalkerGFMC(parentState,NWalkers,NSteps÷nBra,nBra,GuidingWaveFunction,λ)
##
h5open(outfile,"w") do f
    f["TotalWeights",blosc=9] = results.TotalWeights
    f["energies",blosc=9] = results.energies
    f["SaveConfigs",blosc=9] = results.SaveConfigs
    f["reconfigurationTable",blosc=9] = results.reconfTable
    f["nBra"] = nBra
end