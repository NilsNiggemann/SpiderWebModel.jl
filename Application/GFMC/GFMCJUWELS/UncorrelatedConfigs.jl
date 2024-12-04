#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=configs                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/UncorrelatedConfigs/%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/UncorrelatedConfigs.jl ${SLURM_ARRAY_TASK_ID}
exit
=#
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
import SpiderWebModel as SW
using HDF5
using SpiderWebModel.Statistics
i_arg = isinteractive() ? 12 : parse(Int, ARGS[1])

μs = (0.0,0.3,0.4,0.5,0.8,1.0)
# μs = 0.2:0.025:0.45

Ls = (8,12,)
jobs_array = [(;L,μ) for L in Ls for μ in μs]

# μs = μs[1:2:end]
# μs = μs[2:2:end]

(;L,μ) = jobs_array[i_arg]
nBra = 1
τ = 0.2+ 0.1μ
##
NSteps = 20
NConfigs = 8_000
equilibration_steps = 500
pre_equilibration_steps = 0
NWalkers = round(Int,48*100*(L/12)^4)
NStepsEnd = 30
scatter_fraction = 0.0
NBins = 2000
stoch_rec_learning_rate = 3e-3
NWalkers_stochRec = Threads.nthreads() 
equilibration_steps_stochRec = equilibration_steps
report_steps_SR = 2
##
# -- debug params --
if isinteractive()
    # L = 12
    NSteps = 10
    NConfigs = 10
    equilibration_steps = 100
    pre_equilibration_steps = 0
    NWalkers = 12
    stoch_rec_learning_rate = 1e-3
    NStepsEnd = 40
    NBins = 20
    NWalkers_stochRec = Threads.nthreads()
end

##
function get_S_condensate!(S)
    S .= SW.periodicStateDenseLoops(size(S,1))
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
ψGSymm = SW.getNonSymmetric(ψG)
SW.rand!(ψGSymm,1e-3)
SRdir = ENV["MYSCRATCH"]*"/Spiderweb/DataStochRec/L=$(L)/periodic_RK_Full_$(SECTOR_NAME)/$(SW.guidingfunc_name(ψG))/mu=$(μ)/"
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

    SW.get_params(ψG) .= stochReconfRes.params
    w_avg_estimate = -stochReconfRes.E0[end]
elseif isfile(outfileSR) && μ != 1.0
    optim_params_steps = h5read(outfileSR,"params_steps")
    optim_params = selectdim(optim_params_steps,SW.arraydim(optim_params_steps),last(size(optim_params_steps)))
    SW.get_params(ψG) .= optim_params
    w_avg_estimate = -h5read(outfileSR,"E0")[end]
else
    ψG = SW.RKFunction()
    w_avg_estimate = 0.
end

println("stochastic reconf done")
flush(stdout)
sleep(2)

#___________Spin-1_______________________
##

CT = SW.ContinuousTimeMethod(τ,nBra,w_avg_estimate,SW.Hxx_RK(μ))
##

outfileTotal = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_EquilConfigs_2/L=$(L)/mu=$(μ)/Spin1GFMC_Eval_periodic_$(SECTOR_NAME)_mu$(μ)_L$(L)_$(i_arg).h5"
mkpath(dirname(outfileTotal))
# rm(outfileTotal,force=true)
@assert !isfile(outfileTotal) "file already exists!"


SavedConfigs = h5open(outfileTotal,"cw") do file 
    SW.createMMapArray(file,"configs",Int8,(L,L,NConfigs))
end

h5write(outfileTotal,"L",L)
h5write(outfileTotal,"Nwalkers",NWalkers)

initializer = SW.WeightedConfigsInitializers([parentState],[1.])
@info "starting run of $NConfigs" L τ nBra NSteps NWalkers

# outFileDir = ENV["MYSCRATCH"]*"/Spiderweb/temp/DataS1_EquilConfigs_2/L=$(L)/mu=$(μ)/"
# mkpath(outFileDir)

for run in 0:NConfigs
    println(run)
    # outfile = joinpath(outFileDir,"$(run)_($i_arg).h5")
    # rm(outfile,force=true)
    @time results = SW.startManyWalkerGFMC(parentState,CT,NWalkers,NSteps,ψG;equilibration_steps,pre_equilibration_steps=0,initializer)
    empty!(initializer.configs)
    empty!(initializer.weights)
    for walker in 1:NWalkers
        newConf = copy(parentState) 
        newConf .= @view results.SaveConfigs[:,:,walker,end]
        push!(initializer.configs,newConf)
        push!(initializer.weights,1.)
    end
    if run != 0
        SavedConfigs[:,:,run] .= @view results.SaveConfigs[:,:,1,end]
    end
    GC.gc()
    # rm(outfile,force=true)
    flush(stdout)

end
println("GFMC done")
