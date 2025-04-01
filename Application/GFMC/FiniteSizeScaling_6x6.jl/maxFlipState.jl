#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:20274147
#SBATCH --job-name=mF6x6
# SBATCH --job-name=tidyup
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=0-8:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/6x6MC/mf%a.out
#SBATCH --partition=normal
# SBATCH --partition=largemem
#SBATCH --ntasks=1
#SBATCH --mem=230GB
# SBATCH --mem=900GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
# module --force purge
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/FiniteSizeScaling_6x6.jl/maxFlipState.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import ThreadPinning
if !isinteractive()
    ThreadPinning.pinthreads(:cores)
end
import SpiderWebModel as SW
using SpiderWebModel.HDF5
using SpiderWebModel.Statistics
i_arg = isinteractive() ? 2 : parse(Int, ARGS[1])

Ls = (24,30,36)
NRuns = 15
RunBatches = 2
# μs = 0.2:0.025:0.45


jobs_array = [(;L,run) for L in Ls for run in 1:RunBatches:NRuns]

# μs = μs[1:2:end]
# μs = μs[2:2:end]

(;L,run) = jobs_array[i_arg]

tauRange = LinRange(200,100,300)
pre_equilibration_steps = 500_000_000
numRuns=128*30
NSteps = 2000
mu = 0.0
μ = mu
tau=0.3
##
# -- debug params --
if isinteractive()
    # L = 12
    numRuns=100
    RunBatches = 2
    pre_equilibration_steps = 500_000_00
    tauRange = LinRange(50,10,3)
    # ψG = SW.PlaquetteNumberGuidingFunction(0.5)
    NSteps =50
    tau=0.2
end
##
function get_S_condensate!(S)
    S .= 2SW.periodicState6x6Condensate(size(S,1))
    return S
end

SECTOR_NAME  = "6x6Condensate"

parentState = get_S_condensate!(
    SW.stencilConfig(
        zeros(L,L),1,
        boundaryCondition = :periodic
    )
)

##

for run_num in run:min(NRuns, run+RunBatches)
    outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/MaxFlip/$(SECTOR_NAME)/L=$(L)/"
    mkpath(outfileDIR)

    outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_$(run_num).h5")
    if isfile(outfile)
        println("skipping $outfile")
        continue
    end

    @info "starting run $run_num of $NRuns" L outfile
    
    @time maxConf = SW.findMaxFlipConf(parentState;Nwalkers=1,pre_equilibration_steps,NSteps,scatter_fraction=1.,mu=mu,tauRange,numRuns,ψG=SW.PlaquetteNumberGuidingFunction(0.5))

    h5write(outfile,"L",L)
    h5write(outfile,"maxConf",Array(maxConf))
    h5write(outfile,"Nmoves",SW.NPlaquettes(maxConf))
    GC.gc()
    println("GFMC done")
    println("NPlaquettes: ",SW.NPlaquettes(maxConf))
    flush(stdout)
end