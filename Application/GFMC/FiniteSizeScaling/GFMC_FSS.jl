#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8745821
#SBATCH --job-name=smCTRKS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=5-00:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK_FSS/%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/FiniteSizeScaling/GFMC_FSS.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = isinteractive() ? 6 : parse(Int, ARGS[1])

μs = -0.1:0.05:1.1
# μs = 0.2:0.025:0.45

Ls = (20,24,28,32,36)
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
NSteps = 10_000
NBinsEval = 1
NRuns = 14
equilibration_steps = 1200
pre_equilibration_steps = 50_000
NWalkers = round(Int,128*60*(L/28)^2)
scatter_fraction = 0.5
NStepsEnd = 1500
NBins = 1500
stoch_rec_learning_rate = 2e-4
NWalkers_stochRec = Threads.nthreads() * 4
equilibration_steps_stochRec = equilibration_steps
report_steps_SR = 5
##
# -- debug params --
if isinteractive()
    L = 12
    NSteps = 1000
    equilibration_steps = 10
    pre_equilibration_steps = 1000
    NWalkers = 12
    stoch_rec_learning_rate = 1e-4
    NRuns = 2
    NStepsEnd = 30
    NBins = 200
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

ψG = SW.JastrowFunction(parentState)
ψGSymm = SW.getNonSymmetric(ψG)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=($L)/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

outfileSR = if isempty(SRoutfiles)
    joinpath(SRdir,"StochRec_L=$(L)_tau=$(τ)_NW=$(NWalkers)_mu=$(μ).h5")
else
    last(SRoutfiles)
end

CT = SW.ContinuousTimeMethod(τ,1,-length(parentState)*0.266*(1-μ),SW.Hxx_RK(μ))
CT_stochRec = SW.ContinuousTimeMethod(100τ,1,-length(parentState)*0.266*(1-μ),SW.Hxx_RK(μ))

# optimize starting
if !isfile(outfileSR) && μ != 1.0
    @info "starting run" L τ nBra NSteps NWalkers_stochRec outfileSR
    stochReconfRes = SW.stochastic_reconfiguration(parentState,CT_stochRec,NStepsEnd,ψGSymm,NBins,stoch_rec_learning_rate,SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,reconfigure = false,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps=100_000,scatter_fraction,outfile=outfileSR,reset = false,report_steps = report_steps_SR)
end

println("stochastic reconf done")
flush(stdout)

#___________Spin-1_______________________
##

optim_params_steps = h5read(outfileSR,"params_steps")
optim_params = selectdim(optim_params_steps,SW.arraydim(optim_params_steps),last(size(optim_params_steps)))
# optim_params = h5read(outfileSR,"params_steps")[:,:,end]

@assert !iszero(optim_params)
if μ != 1.0
    SW.get_params(ψG) .= optim_params
    # ψG = SW.PlaquetteNumberGuidingFunction(only(unique(optim_params[1])))
    w_avg_estimate = -h5read(outfileSR,"E0")[end]
else
    ψG = SW.RKFunction()
    w_avg_estimate = 0.
end

CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))
##

outfileDIR_init = let
    pth = ENV["MYSCRATCH"]*"/Spiderweb/DataTemp/L=($L)/periodic_RK_Full_$(SECTOR_NAME)/mu=$(μ)/"
    rm(pth,force=true,recursive=true)
    mkpath(pth)

    mktempdir(pth)
end

initializer = getInitializer(parentState,μ,ψG;NWalkers=3*NWalkers,NSteps = 400,OptIndep = 20,outfileDIR=outfileDIR_init)

##
outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/L=$(L)/periodic_RK_Full_$(SECTOR_NAME)/mu_$(μ)/"

Threads.@threads for run in 1:NRuns
    outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(run).h5")
    mkpath(dirname(outfile))

    if !isfile(outfile)
        @info "starting run $run of $NRuns" L τ nBra NSteps NWalkers outfile

        @time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile,initializer)
    end
    flush(stdout)

end
##
println("GFMC done")
flush(stdout)
resFiles = [joinpath(root,file) for (root,_,files) in walkdir(outfileDIR) for file in files]

AllResults = vcat(SW.readResults.(resFiles,NSteps÷NBinsEval)...);
## 
outfileTotal = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_$(SECTOR_NAME)_mu$(μ)_L$(L)_$(i_arg).h5"
mkpath(dirname(outfileTotal))
rm(outfileTotal,force=true)
@assert !isfile(outfileTotal) "file already exists!"
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
    SqsGFMC = SW.getSqsGFMC(AllResults,projectionSteps)
    h5open(outfileTotal,"cw") do file
        file["p_Sq"] = collect(projectionSteps)
        file["SqsGFMC"] = SqsGFMC
    end
    println(L, " ",projectionSteps)
    flush(stdout)
end
rm(outfileDIR,force=true,recursive=true)