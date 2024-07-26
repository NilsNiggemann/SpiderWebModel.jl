#!/bin/bash
#=
#!/bin/bash
# SBATCH --dependency=afterok:8813429
#SBATCH --job-name=CTRKS1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
# SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --time=1-20:00:00
#SBATCH --chdir=/scratch/hpc-prf-pm2frg/niggeni/
#SBATCH --output=/scratch/hpc-prf-pm2frg/niggeni/JobsOutput/Spiderweb/GFMCCTRK/%a_L20.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=220GB
# SBATCH --qos=cont
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
~/.bashrc
julia -O3 -t $SLURM_CPUS_PER_TASK /pc2/groups/hpc-prf-pm2frg/niggeni/Jobs/SpiderWebModel.jl/Application/GFMC/Spin1GFMC_CT_RKL20.jl $SLURM_ARRAY_TASK_ID
exit
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = parse(Int, ARGS[1])

μs = -0.2:0.02:1.2
# μs = μs[1:2:end]
# μs = μs[2:2:end]

μ = μs[i_arg]
L = 20
τ = 0.15
nBra = 1
NSteps = 30_000
NBinsEval = 1
NRuns = 20
equilibration_steps = 10_000
pre_equilibration_steps=0
NWalkers = 128*80
scatter_fraction = 0.0

NStepsstart = 800
NStepsEnd = 1500
NBins = 400
##
# L = 8
# nBra = 1
# NSteps = 1000
# equilibration_steps = 10
# pre_equilibration_steps=100
# NWalkers = 120
# scatter_fraction = 0.6
# NRuns = 1

# NStepsstart = 100
# NBins = 10
##
NWalkers_stochRec = NWalkers ÷ 8
equilibration_steps_stochRec = 100

outfileSR = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataStochRec_periodic_RK_simpl/L=($L)/mu=$(μ)/StochRec_L=$(L)_tau=$(τ)_NW=$(NWalkers)_mu=$(μ)_noscatter.h5"
rm(outfileSR,force=true)

mkpath(dirname(outfileSR))

parentState = SW.stencilConfig(zeros(L,L),1,
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
αstart = 0.13 * (1-μ)
ψG = SW.localPlaquetteGuidingFunction(parentState,αstart)
CT = SW.ContinuousTimeMethod(τ,1,length(parentState)*0.266*(1-μ),SW.Hxx_RK(μ))

##
if !isfile(outfileSR) && μ != 1.0
    @info "starting run" L τ nBra NSteps NWalkers_stochRec outfileSR
    stochReconfRes = SW.stochastic_reconfiguration(parentState,CT,i->min(NStepsstart + 20*i,NStepsEnd),ψG,NBins,i -> 200*min(7,1. +0.02i),SW.IterativeSRSolver();Nwalkers = NWalkers_stochRec,reconfigure = true,rel_tolerance=0.,equilibration_steps=equilibration_steps_stochRec,pre_equilibration_steps=100_000,scatter_fraction,outfile=outfileSR,)
end

println("stochastic reconf done")
flush(stdout)

#___________Spin-1_______________________
##

optim_params = h5read(outfileSR,"params_steps")[:,end]
# optim_params = h5read(infile,"params_steps")[:,:,end]
@assert !iszero(optim_params)
# ψG = SW.FullVariationalGuidingFunction(optim_params)
if μ != 1.0
    ψG = SW.PlaquetteNumberGuidingFunction(only(unique(optim_params[1])))
    w_avg_estimate = -h5read(outfileSR,"E0")[end]
else
    ψG = SW.RKFunction()
    w_avg_estimate = 0.
end

parentState = SW.stencilConfig(zeros(L,L),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))


outfileDIR = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK/L=$(L)/$i_arg/"

for run in 1:NRuns
    outfile = joinpath(outfileDIR,"Spin1GFMC_L=$(L)_tau=$(τ)_NSteps=$(NSteps)_NW=$(NWalkers)_mu=$(μ)_$(run).h5")
    mkpath(dirname(outfile))

    if !isfile(outfile)
        @info "starting run $run of $NRuns" L τ nBra NSteps NWalkers outfile

        @time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps,scatter_fraction,outfile)
    end
    flush(stdout)

end
##
println("GFMC done")
flush(stdout)
resFiles = [joinpath(root,file) for (root,_,files) in walkdir(outfileDIR) for file in files]

AllResults = vcat(SW.readResults.(resFiles,NSteps÷NBinsEval)...);
## 
outfileTotal = "/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataS1_CT_RK/eval/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_mu$(μ)_L$(L)_$(i_arg).h5"
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
    
    Threads.@threads for projectionSteps in (2,50,250,500,750,1000,1250)
    # for projectionSteps in (20,40)
        SqsGFMC = stack(SW.getSqsGFMC(AllResults,projectionSteps),dims=3)
        h5open(outfileTotal,"cw") do file
            file["SqsGFMC/$projectionSteps"] = SqsGFMC
        end
        println(L, " ",projectionSteps)
        flush(stdout)
    
    end
end
