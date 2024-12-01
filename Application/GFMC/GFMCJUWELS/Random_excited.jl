#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=excited                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/Random_excited/%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/Random_excited.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = isinteractive() ? 2 : parse(Int, ARGS[1])


μs = (0.0,0.7,0.8,0.9)
# μs = 0.2:0.025:0.45

Ls = (20,24,28)
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
NRuns = 12
equilibration_steps = 1200
pre_equilibration_steps = 50_000
NWalkers = round(Int,48*60*(L/28)^2)
scatter_fraction = 0.5
NStepsEnd = 30
NBins = 1200
stoch_rec_learning_rate = 1e-3
NWalkers_stochRec = Threads.nthreads()
equilibration_steps_stochRec = equilibration_steps
report_steps_SR = 1
##
# -- debug params --
if isinteractive()
    # L = 12
    NSteps = 500
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
function shuffleSector!(S,N;maxiter = 10N)

    flips = (:diag,:anti)
    # flips = (:row,:col,:diag,:anti)
    spin = SW.getSpin(S)
    function transformSpins!(vec,sgn)
        any(==(sgn*spin),vec) && return
        vec .+= sgn
        return
    end

    function applyMove!(S)
        whichflip = rand(flips)
        sgn = -2sign(sum(S))
        if sgn == 0 
            sgn = rand((-2,2))
        end
        if whichflip == :row
            i = rand(axes(S,1))
            row = @view S[i,begin+iseven(i):2:end]
            transformSpins!(row,sgn)
        elseif whichflip == :col
            j = rand(axes(S,2))
            col = @view S[begin+iseven(j):2:end,j]
            transformSpins!(col,sgn)
        elseif whichflip == :diag
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,1,true)
            transformSpins!(diag,sgn)
        elseif whichflip == :anti
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,-1,true)
            transformSpins!(diag,sgn)
        end
        return 
    end

    for iteration in 1:N
        applyMove!(S)
    end
    for iter in N:maxiter 
        applyMove!(S)
        sum(S) == 0 && break
    end
    # S .= 2S
    SW.fulFillsConstraint(S) || error("Constraint not fulfilled")
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

parentState = 
SW.stencilConfig(
    zeros(L,L),1,
    boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
    )


SW.Random.seed!(i_arg)
NumShuffles = 20000000 + rand(1:20000000)
SECTOR_NAME  = "Random_shuffle_$(NumShuffles)_seed_$(i_arg)"

shuffleSector!(parentState,NumShuffles)

##

# initializer = getInitializer(parentState,μ;NWalkers,NSteps = 1,OptIndep = 2)
##
# rm(outfileSR,force=true)

ψG = SW.SimpleJastrowFunction(parentState)
ψGSymm = SW.getNonSymmetric(ψG)
SW.rand!(ψGSymm,1e-3)

SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$L/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
mkpath(SRdir)
SRoutfiles = readdir(SRdir,join=true)

outfileSR = if isempty(SRoutfiles)
    joinpath(SRdir,"StochRec_L=$(L)_tau=$(τ)_NW=$(NWalkers)_mu=$(μ).h5")
else
    last(SRoutfiles)
end

CT = SW.ContinuousTimeMethod(τ,w_avg_estimate = length(parentState)*0.21*(1-μ),Hxx = SW.Hxx_RK(μ))
CT_stochRec = SW.ContinuousTimeMethod(100τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)
##
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

optim_params = let
    e0 = h5read(outfileSR,"E0") 
    idxzero = findfirst(iszero,e0)
    if isnothing(idxzero)
        optim_params_steps[:,end]
    else
        optim_params_steps[:,idxzero-1]
    end
end
# selectdim(optim_params_steps,SW.arraydim(optim_params_steps),last(size(optim_params_steps)))
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

##
outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/L=$(L)/periodic_RK_Full_$(SECTOR_NAME)/mu_$(μ)/"
outfileTotal = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/eval/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_$(SECTOR_NAME)_mu$(μ)_L$(L)_$(i_arg).h5"
mkpath(dirname(outfileTotal))
# rm(outfileTotal,force=true)
@assert !isfile(outfileTotal) "file already exists!"

initializer = getInitializer(parentState,μ,ψG;NWalkers=3*NWalkers,NSteps = 400,OptIndep = 20,outfileDIR=outfileDIR_init)

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
##
if isinteractive()
    rm(outfileDIR_init,force=true,recursive=true)
    rm(outfileSR,force=true,recursive=true)
    rm(outfileTotal,force=true)
end