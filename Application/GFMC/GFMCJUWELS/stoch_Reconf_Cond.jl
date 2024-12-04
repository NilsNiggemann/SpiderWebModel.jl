#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=SR                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/StochReconf/%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/stoch_Reconf_Cond.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = isinteractive() ? 1 : parse(Int, ARGS[1])


μ = 0.7
# μs = 0.2:0.025:0.45

Ls = (32,36,40)
sector_nums = [2,6,8]
jobs_array = [(L,sector_num) for L in Ls, sector_num in sector_nums]
L,sector_num = jobs_array[i_arg]
nBra = 1
τ = 0.10+ 0.1μ
##
# L = 32
# NSteps = 12_000
# NBinsEval = 1
# NRuns = 14
# equilibration_steps = 800
# pre_equilibration_steps = 50_000
# NWalkers = 128*40
# scatter_fraction = 0.0
# NStepsstart = 1000
# NStepsEnd = 2000
# NBins = 400
##
SW.LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
equilibration_steps = 5200
pre_equilibration_steps = 50_000
scatter_fraction = 0.5
NStepsEnd = 40
NBins = 2000
stoch_rec_learning_rate = 1e-3
NWalkers_stochRec = Threads.nthreads()*4
equilibration_steps_stochRec = equilibration_steps
report_steps_SR = 5
##
# -- debug params --
if isinteractive()
    L = 12
    NSteps = 500
    equilibration_steps = 10
    pre_equilibration_steps = 1000
    NWalkers = 12
    stoch_rec_learning_rate = 1e-6
    NRuns = 2
    NStepsEnd = 10
    NBins = 10
    NWalkers_stochRec = Threads.nthreads() 
end

##

function get_S_condensate!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
    return S
end

SECTOR_NAME  = "4x4_$(sector_num)"

parentState = SW.stencilConfig(
    zeros(L,L),1,
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
parentState .= SW.get4x4PeriodicState(size(parentState,1),sector_num)
##

ψG = SW.SimpleJastrowFunction(parentState)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([-2,2],[2,2]),parentState)
SW.rand!(ψGSymm,1e-3)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$L/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)

getOutfilename(i) = joinpath(SRdir,"StochRec_L=$(L)_tau=$(τ)_NW=$(NWalkers_stochRec)_mu=$(μ)_$(i).h5")

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate = length(parentState)*0.21*(1-μ),Hxx = SW.Hxx_RK(μ))
CT_stochRec = SW.ContinuousTimeMethod(180τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)
##
function findNonZeroParams(filename)
    e0 = h5read(filename,"E0")
    idxzero = findfirst(iszero,e0)
    pars = if isnothing(idxzero)
        h5read(filename,"params_steps")[:,end]
    else
        h5read(filename,"params_steps")[:,idxzero-1]
    end
    @assert !iszero(pars)
    return pars
end
# optimize starting
outfileSR = getOutfilename(i_arg)
while isfile(outfileSR)
    global i_arg += 1
    global outfileSR = getOutfilename(i_arg)
end
if isfile(getOutfilename(i_arg-1))
    SW.get_params(ψG) .= findNonZeroParams(getOutfilename(i_arg-1))
end
@info "starting run" L τ nBra NBins NWalkers_stochRec outfileSR

stochReconfRes = SW.stochastic_reconfiguration(parentState,CT_stochRec,NStepsEnd,ψGSymm,NBins,stoch_rec_learning_rate,SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,reconfigure = false,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps=100_000,scatter_fraction,outfile=outfileSR,reset = false,report_steps = report_steps_SR)

println("stochastic reconf done")
flush(stdout)