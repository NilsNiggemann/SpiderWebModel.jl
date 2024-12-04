#!/bin/bash
#=
#!/bin/bash

#SBATCH --account=pmfrg

#SBATCH --job-name=2ECompGFMC                 # replace name
#SBATCH --export=ALL,JULIA_EXCLUSIVE=1
#SBATCH --mail-user=nils.niggemann@fu-berlin.de  # replace email address
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90GB         # memory , more means less gc time
#SBATCH --time=0-24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=END
#SBATCH --output=/p/project/pmfrg/niggemann1/JobsOutput/Spiderweb/GFMC/EnergySweep/4x4Spin1_%a.out    # File to which standard Out- will be written

jutil env activate -p pmfrg
cd $PROJECT/niggemann1
module --force purge
module load Stages/2024  
module load GCCcore/.12.3.0

module load Julia/1.9.3
export JULIA_DEPOT_PATH=/p/scratch/pmfrg/niggemann1/.julia/

julia -O3 -t $SLURM_CPUS_PER_TASK /p/project/pmfrg/niggemann1/Jobs/SpiderWebModel.jl/Application/GFMC/GFMCJUWELS/4x4Orders_distr.jl ${SLURM_ARRAY_TASK_ID}
exit
=#

import Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__)
import Pkg
Pkg.instantiate()
Pkg.precompile()
import SpiderWebModel as SW
using SpiderWebModel.HDF5
i_arg = parse(Int, ARGS[1])

function upscale(Conf,L)
    S = SW.stencilConfig(zeros(L,L),SW.getSpin(Conf);
    boundaryCondition = :periodic
    )

    per = SW.PeriodicMatrix(Conf,L,L)
    for I in CartesianIndices(S)
        S[I] = per[I]
    end
    S
end
function findFirstMiIndex(arr)
    minval = arr[begin]
    i_min = firstindex(arr)
    for (i,x) in enumerate(arr)
        if x < minval
            i_min = i
            minval = x
        end
    end
    return i_min
end

function optimizeWF(S,CT)
    CTSR = SW.ContinuousTimeMethod(100*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

    psi = SW.SimpleJastrowFunction(S)

    ψGSymm = SW.symmetrize(psi,SW.TranslationalSymmetry([-2,2],[2,2]),S)

    stochReconfResSymm = SW.stochastic_reconfiguration(S,CTSR,20 ,ψGSymm,1500,8e-3,SW.IterativeSRSolver();Nwalkers = 2*Threads.nthreads(),reconfigure=false,rel_tolerance=0,equilibration_steps=100,pre_equilibration_steps=40_000,
    report_steps = 100,
    reset = false,
    # outfile = "tempSR/SR2.h5"
    )
    SW.get_params(psi) .= stochReconfResSymm.params
    return psi
end

function findEnergies(Configs,CT;Nwalkers = 28*2,NSteps = 2000,equilibration_steps = NSteps ÷6,outfileDIR = nothing)
    en = zeros(length(Configs))
    Δen = zeros(length(Configs))
    getOutfile(i,j) = isnothing(outfileDIR) ? nothing : joinpath(outfileDIR,"$(i)_$(j).h5")

    for i in eachindex(Configs,en)
        S = Configs[i]
        # display(SW.plotApplPlaquettes(S))
        if length(SW.getApplicablePlaquettes(S)) == 0
            en[i] = 0
            Δen[i] = 0
            continue
        end

        ψG = optimizeWF(S,CT)

        NwalkersNew = max(48*2, round(Int,Nwalkers * (sum(SW.getNPlaq(S))/ 48 / 10)^2)*48)
        @info "starting $i" NwalkersNew
        
        results = [SW.startManyWalkerGFMC(S,CT,NwalkersNew,NSteps,ψG;equilibration_steps,pre_equilibration_steps=NSteps,scatter_fraction = 0.8,outfile = getOutfile(i,i_st)) for i_st in 1:6]
        GC.gc()
        energies = SW.getEnergies.(results,1,min(NSteps÷3, 300))
        energiesMean = SW.mean.(energies)
        energiesStd = SW.std.(energies)
        e0_index = findFirstMiIndex(energiesMean)
        en[i] = energiesMean[e0_index]
        Δen[i] = energiesStd[e0_index]
    end
    return (;en,Δen)
end
function makeConf(UC,L,Spin)
    S = SW.stencilConfig(zeros(L,L),Spin;
    boundaryCondition = :periodic
    )
    S .= SW.getPeriodicState(UC,L,L)
    return S
end  
##
L = 16
Nwalkers = 48*15
NSteps = 6000

reducedConfigs = makeConf.(collect.(eachslice(SW.h5read("../../Data/reducedConfigs.h5","reducedConfigs"),dims=3)),L,1)


mus_sectors = LinRange(-0.15,0.99,30)
mu = mus_sectors[i_arg]
##
CT = SW.ContinuousTimeMethod(0.2,1,-0.266length(reducedConfigs[1]),SW.Hxx_RK(mu))
outfileDIR = ENV["MYSCRATCH"]*"/Spiderweb/DataS1_CT_RK_equil/SectorComp/runs/L=$(L)/mu=$(mu)_$(i_arg)/"
rm(outfileDIR,recursive=true,force=true)
mkpath(outfileDIR)

res = findEnergies(upscale.(reducedConfigs,L),CT;Nwalkers,NSteps,outfileDIR)
##
AllresEn = res.en
AllresΔEn = res.Δen

outfileTotal = "/p/scratch/pmfrg/niggemann1/Spiderweb/DataS1_CT_RK_equil/SectorComp3/L=$(L)/mu=$(mu).h5"
mkpath(dirname(outfileTotal))

h5write(outfileTotal,"energies",AllresEn)
h5write(outfileTotal,"Δenergies",AllresΔEn)
h5write(outfileTotal,"mu",mu)
##
rm(outfileDIR,recursive=true)