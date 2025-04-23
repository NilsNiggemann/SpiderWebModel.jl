#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:21839386
#SBATCH --job-name=LmuStair
# SBATCH --job-name=tidyup
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=0-18:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK_Staircase/%a.out
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

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/StairCase_Sweep.jl/GFMC_Stair_Sq_fine.jl $SLURM_ARRAY_TASK_ID
exit
=#
#29,30,31,32,33,34,35,36,37,38,39,40,41,42,99,100,101,102,103,104,105,106,107,108,109,110,111,112,169,170,171,172,173,174,175,176,177,178,179,180,181,182
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
i_arg = isinteractive() ? 30 : parse(Int, ARGS[1])

# μs = range(0.81,0.89,length = 5)
μs = (0.78,0.81,0.85)
NRuns = 10
# RunBatches = 1
function RunBatchesFunc(L)
    if L == 28
        return 14
    elseif L == 32
        return 7
    else
        return 1
    end
end
# μs = 0.2:0.025:0.45

Ls = (28,32,36)
jobs_array = [(;L,mu,run) for L in Ls for mu in μs for run in 1:RunBatchesFunc(L):NRuns]

# μs = μs[1:2:end]
# μs = μs[2:2:end]
# 1,2,3,6,7,8,9,10,11,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57
(;L,mu,run) = jobs_array[i_arg]
RunBatches = RunBatchesFunc(L)
μ = mu
τ = 0.1+ 0.1μ
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
pre_equilibration_steps = round(Int,500_000_000*(L/36))
NWalkers = round(Int,128*30*(L/24)^4)
NWalkers = (NWalkers - NWalkers%128)
scatter_fraction = 0.7
projection_order = 120
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
SECTOR_NAME  = "StairCase"

parentState = get_S_stair!(
    SW.stencilConfig(
        zeros(L,L),1,
        boundaryCondition=:periodic
    )
)

##

# initializer = getInitializer(parentState,μ;NWalkers,NSteps = 1,OptIndep = 2)
##
# rm(outfileSR,force=true)

ψG = SW.SimpleJastrowFunction(parentState)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$L/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate = length(parentState)*0.21*(1-μ),Hxx = SW.Hxx_RK(μ))
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
    crit1 = (e0diff[end] < 1e-3) 
    crit2 = (ΔEdiff[end] < 5e-4)
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
if μ != 1
    _SR_iteration = 1
    isempty(SRoutfiles) && error("SRoutfiles empty")
    outfileSR = last(SRoutfiles)
    if !convergence_heuristic(outfileSR)
        error("no convergence of SR!")
    end

    idx = findNonZeroEn(outfileSR)
    SW.get_params(ψG) .= findNonZeroParams(outfileSR,idx)
    w_avg_estimate = h5read(outfileSR,"E0")[idx]
else
    w_avg_estimate = 0.0
    ψG = SW.RKFunction()
end
#___________Spin-1_______________________
##

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate,SW.Hxx_RK(μ))

##

##
# initializer = getInitializer(parentState,μ,ψG;NWalkers=NWalkers,NSteps = 100,OptIndep = 6,outfileDIR=outfileDIR_init)
GC.gc()

initializer = SW.CombinedInitializer(
    SW.UnguidedWalkInitializer(pre_equilibration_steps,scatter_fraction), 
    SW.StochasticResettingInitializer(LinRange(200,CT.τ,60),CT,800.,parentState)
)
for run_num in run:min(NRuns, run+RunBatches)
    outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/$(SECTOR_NAME)_fine_2/L=$(L)/mu=$(μ)/$run_num/"
    mkpath(outfileDIR)

    outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(run_num).h5")
    if isfile(outfile)
        println("skipping $outfile")
        continue
    end

    @info "starting run $run_num of $NRuns" L τ NSteps NWalkers outfile
    h5write(outfile,"L",L)
    h5write(outfile,"mu",μ)

    @time results = SW.measure_Sq_GFMC(parentState,CT,NWalkers,NSteps,projection_order,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile,estimate_w_avg = true,initializer)
    # initializer
    GC.gc()
    println("GFMC done")

    flush(stdout)
end