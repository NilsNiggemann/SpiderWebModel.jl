#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8745821
#SBATCH --job-name=smCTRKS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-00:00:00
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
i_arg = parse(Int, ARGS[1])

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
NStepsstart = 1000
NStepsEnd = 1500
NBins = 300

## -- debug params --
# L = 32
# nBra = 1
# NSteps = 1000
# equilibration_steps = 10
# pre_equilibration_steps=100
# NWalkers = 12
# scatter_fraction = 0.6
# NRuns = 1
# NStepsstart = 200
# NStepsEnd = 500
# NBins = 10

# NStepsstart = 100
# NBins = 10
##
NWalkers_stochRec = NWalkers ÷ 4
equilibration_steps_stochRec = 100
##
function get_S_condensate!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
    return S
end

function get_S_stair!(S)
    S .= 2SW.getStaircase(size(S,1))
    return S
end

function getInitializer(S,mu;NWalkers=128,NSteps = 100,tau =0.3 + 0.6mu,OptIndep = 10,outfileDIR=nothing)
    μ = mu
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
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

outfileDIR_init = let
    pth = ENV["MYSCRATCH"]*"/Spiderweb/DataTemp/L=($L)/periodic_RK_Full_$(SECTOR_NAME)/mu=$(μ)/"
    mkpath(pth)

    mktempdir(pth)
end

initializer = getInitializer(parentState,μ;NWalkers=3*NWalkers,NSteps = 400,OptIndep = 20,outfileDIR=outfileDIR_init)
# initializer = getInitializer(parentState,μ;NWalkers,NSteps = 1,OptIndep = 2)
##
# rm(outfileSR,force=true)
SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=($L)/periodic_RK_Full_$(SECTOR_NAME)/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

outfileSR = if isempty(SRoutfiles)
    joinpath(SRdir,"StochRec_L=$(L)_tau=$(τ)_NW=$(NWalkers)_mu=$(μ).h5")
else
    last(SRoutfiles)
end

αstart = 0.15*(1-μ)
ψG = SW.orderGuidingFunction(parentState,0.12)
CT = SW.ContinuousTimeMethod(τ,1,-length(parentState)*0.266*(1-μ),SW.Hxx_RK(μ))
ψGSymm = SW.symmetrize(parentState,ψG,(4,4))

##

# optimize starting
if !isfile(outfileSR) && μ != 1.0
    @info "starting run" L τ nBra NSteps NWalkers_stochRec outfileSR
    stochReconfRes = SW.stochastic_reconfiguration(parentState,CT,i->min(NStepsstart + 20*i,NStepsEnd),ψGSymm,NBins,i -> 1*max(0.3,0.6 -0.002i),SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,reconfigure = true,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps=100_000,scatter_fraction,outfile=outfileSR,initializer,reset = false)
end

println("stochastic reconf done")
flush(stdout)

#___________Spin-1_______________________
##

optim_params_steps = h5read(outfileSR,"params_steps")
optim_params = selectdim(optim_params_steps,SW.arraydim(optim_params_steps),last(size(optim_params_steps)))
# optim_params = h5read(outfileSR,"params_steps")[:,:,end]
Nplaquettes = length(collect(SW.plaquetteIterator(parentState)))

@assert !iszero(optim_params)
if μ != 1.0
    ψG.params .= optim_params
    # ψG = SW.PlaquetteNumberGuidingFunction(only(unique(optim_params[1])))
    w_avg_estimate = -h5read(outfileSR,"E0")[end]
else
    ψG = SW.RKFunction()
    w_avg_estimate = 0.
end

CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))

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
    
    projectionSteps = 1:5:150
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