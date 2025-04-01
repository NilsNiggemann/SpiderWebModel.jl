#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:16952754
#SBATCH --job-name=SR_gen
# SBATCH --job-name=tidyup
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-20:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/SR/%a.out
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

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/GenericSector.jl/StochRec.jl $SLURM_ARRAY_TASK_ID
exit
=#
i_arg = isinteractive() ? 2 : parse(Int, ARGS[1])
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())


import ThreadPinning
if !isinteractive()
    ThreadPinning.pinthreads(:cores)
end
import SpiderWebModel as SW
using SpiderWebModel.HDF5
using SpiderWebModel.Statistics

μs = (0.6,0.7,0.8,0.9,0.95)
confs = h5read(ENV["MYSCRATCH"]*"confs/confs_fixed_conserved_2.h5","confs")
L = size(confs)[1]
numConfs = size(confs)[3]
jobs_array = [(;confnum,mu) for mu in μs for confnum in 1:5]

# μs = μs[1:2:end]
# μs = μs[2:2:end]

(;confnum,mu) = jobs_array[i_arg]
μ = mu
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
pre_equilibration_steps = 1_000_000_000
scatter_fraction = 1.0
NStepsEnd = 40
NBins = 2000
stoch_rec_learning_rate = 2e-2

NWalkers_stochRec = Threads.nthreads() * 2
equilibration_steps_stochRec = 500
report_steps_SR = 5
##
# -- debug params --
if isinteractive()
    # L = 12
    equilibration_steps = 10
    pre_equilibration_steps = 1000
    NWalkers = 12
    NBins = 10

    stoch_rec_learning_rate = 1e-6
    NStepsEnd = 10
    NWalkers_stochRec = Threads.nthreads() ÷ 2
end

##
SECTOR_NAME  = "RandConf_$(confnum)"

parentState = SW.stencilConfig(
        zeros(L,L),1,
        boundaryCondition = :periodic
)
parentState .= confs[:,:,confnum]
##

ψG = SW.SimpleJastrowFunction(parentState)
Symmetry = SW.TranslationalSymmetry(SW.SA[1,-1],SW.SA[1,1])
# Symmetry = SW.SymmetryGroup(SW.ExchangeSymmetry())
# ψGSymm = SW.getNonSymmetric(ψG)
ψGSymm = SW.symmetrize(ψG,Symmetry,parentState)
# ψGSymm = SW.getNonSymmetric(ψG)
SW.rand!(ψGSymm,1e-5)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$L/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

getOutfilename(i) = joinpath(SRdir,"StochRec_L=$(L)_tau=$(CT_stochRec.τ)_NW=$(NWalkers_stochRec)_mu=$(μ)_$(i).h5")

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate = length(parentState)*0.21*(1-μ),Hxx = SW.Hxx_RK(μ))
CT_stochRec = SW.ContinuousTimeMethod(10τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)
##
function findNonZeroEn(filename)
    e0 = h5read(filename,"E0")
    idxzero = findfirst(iszero,e0)
    if isnothing(idxzero)
        return lastindex(e0)
    else
        return idxzero-1
    end
end

function findNonZeroParams(filename,idxzero = findNonZeroEn(filename))
    pars = h5read(filename,"params_steps")[:,idxzero]
    @assert !iszero(pars)
    return pars
end

function movingaverage(X::AbstractVector,numofele::Int)
    BackDelta = div(numofele,2) 
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = lastindex(X)
    Y = similar(X)
    for n = eachindex(X)
        lo = max(firstindex(X),n - BackDelta)
        hi = min(len,n + ForwardDelta)
        @views Y[n] = mean(X[lo:hi])
    end
    return Y
end

function convergence_heuristic(filename)
    idx = findNonZeroEn(filename)
    e0 = movingaverage(h5read(filename,"E0")[begin:idx],300)
    ΔE = movingaverage(h5read(filename,"ΔE")[begin:idx],300)

    e0diff = abs.(diff(e0))
    ΔEdiff = abs.(diff(ΔE))
    crit1 = (e0diff[end] < 1e-4) 
    crit2 = (ΔEdiff[end] < 1e-4)
    return crit1 && crit2
end
# optimize starting
_SR_iteration = 1
outfileSR = getOutfilename(_SR_iteration)
if !isempty(SRoutfiles)
    
    outfileSR = last(SRoutfiles)
    idx = findNonZeroEn(outfileSR)
    SW.get_params(ψG) .= findNonZeroParams(outfileSR,idx)
    SW.enforceSymmetries!(ψGSymm)

    iternum = length(SRoutfiles)
    if !convergence_heuristic(outfileSR) || true
        global _SR_iteration = iternum + 1
        global outfileSR = getOutfilename(_SR_iteration)
    end
end
##
# optimize starting


if !isfile(outfileSR)
    @info "starting run" L τ NStepsEnd NWalkers_stochRec outfileSR
    stochReconfRes = SW.stochastic_reconfiguration(parentState,CT_stochRec,NStepsEnd,ψGSymm,NBins,stoch_rec_learning_rate,SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps,scatter_fraction,outfile=outfileSR,reset = false,report_steps = report_steps_SR)

    # ψG = SW.PlaquetteNumberGuidingFunction(only(unique(optim_params[1])))
end

println("stochastic reconf done")
println("convergence: ", convergence_heuristic(outfileSR))
flush(stdout)
