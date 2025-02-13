#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:16952754
#SBATCH --job-name=L40CTRKS1
# SBATCH --job-name=tidyup
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-20:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK_FSS/%a.out
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

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/FiniteSizeScaling_directSq.jl/GFMC_FSS_Sq.jl $SLURM_ARRAY_TASK_ID
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
i_arg = isinteractive() ? 80 : parse(Int, ARGS[1])

μs = -0.1:0.05:1.0
NRuns = 14
RunBatches = 5
# μs = 0.2:0.025:0.45

Ls = (32,36,40)
jobs_array = [(;L,μ,run) for L in Ls for μ in μs for run in 1:RunBatches:NRuns]

# μs = μs[1:2:end]
# μs = μs[2:2:end]

(;L,μ,run) = jobs_array[i_arg]
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
NSteps = 10_000
NBinsEval = 1
equilibration_steps = 1000
pre_equilibration_steps = 50_000
NWalkers = round(Int,128*20*(L/24)^4)
if 0.1<= μ <= 0.5
    NWalkers *= 4
end
NWalkers = (NWalkers - NWalkers%128)
scatter_fraction = 0.5
projection_order = 150
##
# -- debug params --
if isinteractive()
    # L = 12
    NSteps = 100
    equilibration_steps = 10
    pre_equilibration_steps = 1000
    NWalkers = 12
    NRuns = 2
    projection_order = 20
    NStepsEnd = 10
end

##
function get_S_condensate!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
    return S
end

function get_S_stair!(S)
    S .= 4SW.getStairCase(size(S,1))
    return S
end
##
SECTOR_NAME  = "Condensate"

parentState = get_S_condensate!(
    SW.stencilConfig(
        zeros(L,L),1,
        boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
)

##

# initializer = getInitializer(parentState,μ;NWalkers,NSteps = 1,OptIndep = 2)
##
# rm(outfileSR,force=true)

ψG = SW.SimpleJastrowFunction(parentState)
Symmetry = SW.TranslationalSymmetry(SW.SA[2,2],SW.SA[-2,2])
# Symmetry = SW.SymmetryGroup(SW.ExchangeSymmetry())
# ψGSymm = SW.getNonSymmetric(ψG)
ψGSymm = SW.symmetrize(ψG,Symmetry,parentState)
# ψGSymm = SW.getNonSymmetric(ψG)
SW.rand!(ψGSymm,1e-3)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$L/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

getOutfilename(i) = joinpath(SRdir,"StochRec_L=$(L)_tau=$(CT_stochRec.τ)_NW=$(NWalkers_stochRec)_mu=$(μ)_$(i).h5")

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate = length(parentState)*0.21*(1-μ),Hxx = SW.Hxx_RK(μ))
CT_stochRec = SW.ContinuousTimeMethod(150τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)
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
    # @info "" e0diff[end] ΔEdiff[end]
    crit1 = (e0diff[end] < 1e-3) 
    crit2 = (ΔEdiff[end] < 1e-3)
    return crit1 && crit2
    # fig = Figure()
    # ax = Axis(fig[1,1],title = "$crit1 , $crit2")
    # lines!(ax,e0diff)
    # lines!(ax,ΔEdiff)
    # display(fig)
    # display(plotVarEn(filename,plotDiff=false))
    # return fig

end
# optimize starting
_SR_iteration = 1
isempty(SRoutfiles) && error("SRoutfiles not empty")
outfileSR = last(SRoutfiles)
if !convergence_heuristic(outfileSR)
    @warn "no full convergence of SR!"
end

idx = findNonZeroEn(outfileSR)
SW.get_params(ψG) .= findNonZeroParams(outfileSR,idx)
w_avg_estimate = h5read(outfileSR,"E0")[idx]
#___________Spin-1_______________________
##

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate,SW.Hxx_RK(μ))

##

##
# initializer = getInitializer(parentState,μ,ψG;NWalkers=NWalkers,NSteps = 100,OptIndep = 6,outfileDIR=outfileDIR_init)
GC.gc()


for run_num in run:min(NRuns, run+RunBatches)
    outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/L=$(L)/mu=$(μ)/$run_num/"
    mkpath(outfileDIR)

    outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(run_num).h5")
    if isfile(outfile)
        println("skipping $outfile")
        continue
    end

    @info "starting run $run_num of $NRuns" L τ NSteps NWalkers outfile
    h5write(outfile,"L",L)
    h5write(outfile,"mu",μ)

    @time results = SW.measure_Sq_GFMC(parentState,CT,NWalkers,NSteps,projection_order,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile,estimate_w_avg = true)
    # initializer
    GC.gc()
    println("GFMC done")

    flush(stdout)
end