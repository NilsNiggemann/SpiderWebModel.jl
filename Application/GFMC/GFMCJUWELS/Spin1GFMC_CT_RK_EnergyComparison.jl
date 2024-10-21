#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=ECompGFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/EnergySweep/Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/Spin1GFMC_CT_RK_EnergyComparison.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = parse(Int, ARGS[1])

μs = -0.3:0.1:1.1
Sectors = [
    "S0",
    "DenseLoops",
    "Diag",
]
# μs = μs[1:2:end]
# μs = μs[2:2:end]

jobarray = [(;μ,SECTOR_NAME) for μ in μs for SECTOR_NAME in Sectors]

(;μ,SECTOR_NAME) = jobarray[i_arg]
nBra = 1
τ = 0.10
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
L = 24
NSteps = 6_000
NBinsEval = 1
NRuns = 12
equilibration_steps = 800
pre_equilibration_steps = 30_000
NWalkers = 128*10
scatter_fraction = 0.0
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
function get_S_DenseLoops!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
    return S
end

function getInitializer(S,mu;NWalkers=128,NSteps = 100,tau = 1. + 2mu,OptIndep = 10)
    μ = mu
    ψG = SW.PlaquetteNumberGuidingFunction(0.15*(1-μ))
    CTFindOpt = SW.ContinuousTimeMethod(tau,1,(1-μ)* 0.266*length(S),SW.Hxx_RK(μ))
            
    @time OptimStart = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CTFindOpt,NWalkers,NSteps,ψG;equilibration_steps=1,pre_equilibration_steps=60_000,scatter_fraction=0.5) for _ in 1:OptIndep])
    initializer = SW.WeightedConfigsInitializers(OptimStart)
    return initializer
end
##

function get_S_stair(S)
    S_staircase = copy(S)
    S_staircase .= 4*SW.getStairCase(size(S,1))
    return S_staircase
end

parentState = let
    
    S0 = SW.stencilConfig(
        zeros(L,L),1,
        boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )
    # "S0",
    # "StairCase",
    # "DenseLoops",
    # "Diag"
    if SECTOR_NAME == "S0"
    elseif SECTOR_NAME == "DenseLoops"
        get_S_DenseLoops!(S0)
    elseif SECTOR_NAME == "Diag"
        S0 .= SW.periodicStateDiag(size(S0,1))
    end
    S0
end

##

initializer = getInitializer(parentState,μ;NWalkers,NSteps = 100,OptIndep = 20)
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

for run in 1:NRuns
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
@assert !isfile(outfileTotal) "file already exists!"
function getEns(results)
    en = [SW.getEnergies(res.TotalWeights,res.energies,1,1000÷res.nBra) for res in results]
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
    end
    
    Threads.@threads for projectionSteps in (1,10,20,50,100,250,300)
    # for projectionSteps in (20,40)
        SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
        h5open(outfileTotal,"cw") do file
            file["SqsGFMC/$projectionSteps"] = SqsGFMC
        end
        println(L, " ",projectionSteps)
        flush(stdout)
    
    end
end
