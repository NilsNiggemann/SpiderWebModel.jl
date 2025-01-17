#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:16952754
# SBATCH --job-name=L28CTRKS1
#SBATCH --job-name=tidyup
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-20:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK_FSS/%a.out
# SBATCH --partition=normal
#SBATCH --partition=largemem
#SBATCH --ntasks=1
# SBATCH --mem=230GB
#SBATCH --mem=900GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
# module --force purge
module load lang/JuliaHPC/1.10.1-foss-2022a-CUDA-11.7.0

julia -O3 -t $SLURM_CPUS_PER_TASK --heap-size-hint=210G /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/FiniteSizeScaling/GFMC_FSS.jl $SLURM_ARRAY_TASK_ID
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
i_arg = isinteractive() ? 84 : parse(Int, ARGS[1])

μs = -0.1:0.05:1.1
# μs = 0.2:0.025:0.45

Ls = (16,20,24,28,32,36)
jobs_array = [(;L,μ) for L in Ls for μ in μs]

# μs = μs[1:2:end]
# μs = μs[2:2:end]

(;L,μ) = jobs_array[i_arg]
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
NSteps = 5_000
NBinsEval = 1
NRuns = 14
equilibration_steps = 1000
pre_equilibration_steps = 50_000
NWalkers = round(Int,128*50*(L/24)^4)
if 0.1<= μ <= 0.5
    NWalkers *= 4
end
NWalkers = (NWalkers - NWalkers%128)
scatter_fraction = 0.5
NStepsEnd = 40
NBins = 4000
stoch_rec_learning_rate = 1e-2
if L == 28
    stoch_rec_learning_rate = 6e-3
end
NWalkers_stochRec = Threads.nthreads() * 3
equilibration_steps_stochRec = equilibration_steps
report_steps_SR = 1
##
# -- debug params --
if isinteractive()
    # L = 12
    NSteps = 1000
    equilibration_steps = 10
    pre_equilibration_steps = 1000
    NWalkers = 12
    stoch_rec_learning_rate = 1e-6
    NRuns = 2
    NStepsEnd = 10
    NBins = 10
    NWalkers_stochRec = Threads.nthreads() ÷ 2
end

##
function get_S_condensate!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
    return S
end

function get_S_stair!(S)
    S .= 2SW.getStaircase(size(S,1))
    return S
end

function getInitializer(S,mu,ψG;NWalkers=128,NSteps = 100,tau =0.5 + 0.4mu,OptIndep = 10,outfileDIR=nothing)
    μ = mu
    CTFindOpt = SW.ContinuousTimeMethod(tau,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))
    
    getOutf(i) = isnothing(outfileDIR) ? nothing : joinpath(outfileDIR,"$i.h5")

    @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NWalkers,NSteps,ψG;equilibration_steps=1,pre_equilibration_steps=60_000,scatter_fraction=0.5,outfile=getOutf(i)) for i in 1:OptIndep])
    initializer = SW.WeightedConfigsInitializers(OptimStart)
    return initializer
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
_SR_iteration = 1
outfileSR = getOutfilename(_SR_iteration)
if !isempty(SRoutfiles)
    
    outfileSR = last(SRoutfiles)
    if !convergence_heuristic(outfileSR)
        global _SR_iteration += 1
        global outfileSR = getOutfilename(_SR_iteration)
    end
end

if isfile(getOutfilename(_SR_iteration-1))
    idx = findNonZeroEn(filename)
    SW.get_params(ψG) .= findNonZeroParams(getOutfilename(_SR_iteration-1),idx)
    SW.enforceSymmetry!(ψGSymm)
end
##
# optimize starting
if μ==1.0
    ψG = SW.RKFunction()
    w_avg_estimate = 0.
elseif !isfile(outfileSR)
    @info "starting run" L τ nBra NSteps NWalkers_stochRec outfileSR
    stochReconfRes = SW.stochastic_reconfiguration(parentState,CT_stochRec,NStepsEnd,ψGSymm,NBins,stoch_rec_learning_rate,SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,reconfigure = false,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps=100_000,scatter_fraction,outfile=outfileSR,reset = false,report_steps = report_steps_SR)

        
    SW.get_params(ψG) .= stochReconfRes.params
    # ψG = SW.PlaquetteNumberGuidingFunction(only(unique(optim_params[1])))
    w_avg_estimate = -stochReconfRes.E0[end]

else
    idx = findNonZeroEn(outfileSR)
    SW.get_params(ψG) .= findNonZeroParams(outfileSR,idx)
    w_avg_estimate = -h5read(outfileSR,"E0")[idx]
end

println("stochastic reconf done")
println("convergence: ", convergence_heuristic(outfileSR))
flush(stdout)

#___________Spin-1_______________________
##

CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))
##

outfileDIR_init = let
    pth = ENV["MYSCRATCH"]*"/Spiderweb/DataTemp/L=($L)/periodic_RK_Full_$(SECTOR_NAME)/mu=$(μ)/"
    rm(pth,force=true,recursive=true)
    mkpath(pth)

    mktempdir(pth)
end

##

function is_valid_file(filename)
    !isfile(filename) && return false
    allkeys = ["energies","nBra","L","mu","tau","NWalkers","NSteps","SqsGFMC","p_Sq"]
    h5open(filename,"r") do file
        return all(k->haskey(file,k),allkeys)
    end
end

outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/L=$(L)/periodic_RK_Full_$(SECTOR_NAME)/mu_$(μ)/"
outfileTotal = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_$(SECTOR_NAME)_mu$(μ)_L$(L)_$(i_arg).h5"
mkpath(dirname(outfileTotal))
# rm(outfileTotal,force=true)
if !is_valid_file(outfileTotal)
    println("Deleting invalid output file:\n $outfileTotal")
    rm(outfileTotal,force=true)
else 
    println("output $outfileTotal already exists! Exiting...")
    exit(0) 
end

@assert !isfile(outfileTotal) "file already exists!"
##
# initializer = getInitializer(parentState,μ,ψG;NWalkers=NWalkers,NSteps = 100,OptIndep = 6,outfileDIR=outfileDIR_init)
GC.gc()
let 
    Threads.@threads for run in 1:NRuns
        outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(run).h5")
        mkpath(dirname(outfile))

        if !isfile(outfile)
            @info "starting run $run of $NRuns" L τ nBra NSteps NWalkers outfile

            @time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile
            # initializer
            )
        end
        GC.gc()
        flush(stdout)
    end

end

##
println("GFMC done")
flush(stdout)
GC.gc()
resFiles = [joinpath(root,file) for (root,_,files) in walkdir(outfileDIR) for file in files]

AllResults = vcat(SW.readResults.(resFiles,NSteps÷NBinsEval)...);
## 

function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,500÷res.nBra) for res in results]
end
##

let 
    en  = stack(getEns(AllResults))
    L = size(AllResults[1].SaveConfigs,1)
    h5open(outfileTotal,"cw") do file
        file["energies"] = en
        file["nBra"] = AllResults[1].nBra
        file["L"] = L
        file["mu"] = μ
        file["tau"] = τ
        file["NWalkers"] = NWalkers
        file["NSteps"] = NSteps
    end
    
    projectionSteps = 1:150
    # for projectionSteps in (20,40)
    SqsGFMC = SW.getSqsGFMC(AllResults,projectionSteps,useBuffer=true)
    h5open(outfileTotal,"cw") do file
        file["p_Sq"] = collect(projectionSteps)
        file["SqsGFMC"] = SqsGFMC
    end
    println(L, " ",projectionSteps)
    flush(stdout)
end



if is_valid_file(outfileTotal)
    println("Final output created. Deleting temp dir ",outfileDIR)
    rm(outfileDIR,force=true,recursive=true)
else
    println("Final output not created. Keeping temp dir. ",outfileDIR)
end

# rm(outfileDIR,force=true,recursive=true)