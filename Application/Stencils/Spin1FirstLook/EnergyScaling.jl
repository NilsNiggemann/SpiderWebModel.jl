using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

using ThreadPinning
ThreadPinning.pinthreads(:cores)
import Pkg
cd(@__DIR__)
Pkg.activate("../../")
if !Pkg.Operations.is_instantiated(Pkg.Types.EnvCache(Base.active_project()))
    @info "instantiating environment"
    Pkg.instantiate()
end
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using SpiderWebModel.HDF5
include("../plottingUtils.jl")

##
#___________Periodic Boundaries_______________________

function initializeGWF(S,Symmetry) 
    ψG = SW.SimpleJastrowFunction(S)
    ψGSymm = SW.symmetrize(ψG,Symmetry,S)
    SW.rand!(ψGSymm,1e-5)
    return ψGSymm
end

function _optimizeWF(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,NConfs,ψG,NSteps,dt,SW.IterativeSRSolver();Nwalkers = 28*4,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=25,kwargs...)
    if showplot
        display(plotVarEn(stochReconfRes))
    end
    return stochReconfRes
end

function optimizeWF!(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = _optimizeWF(ψG,S,CTSR,nThermal;dt = dt,NSteps = NSteps,NConfs = NConfs,showplot=showplot,kwargs...)

    E0 = stochReconfRes.E0
    idx = 0
    NSteps_base = NSteps
    while E0[end] > E0[begin]
        @warn "Optimization failed, reducing dt and increasing NSteps"
        idx += 1
        # idx > 5 && error("Optimization failed")
        idx > 4 && break
        dt *= 0.75
        NSteps += 2NSteps_base
        SW.rand!(ψG,1e-6)
        stochReconfRes = _optimizeWF(ψG,S,CTSR,nThermal;dt = dt,NSteps = NSteps,NConfs = NConfs,showplot=showplot,kwargs...)
        E0 = stochReconfRes.E0
    end
    SW.get_params(ψG) .= stochReconfRes.params

    return stochReconfRes.E0[end]
end
##

function getEnergy(S,mu,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,dt = 3e-3,NConfs = 20,nopt=1000,pre_init_steps=800)
    ψG = initializeGWF(S,Symmetry)
    SW.rand!(ψG,1e-6)
    tau = 0.1 + 0.1*mu

    CTSR = SW.ContinuousTimeMethod(30*tau,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(mu))

    @time E0_est = optimizeWF!(ψG,S,CTSR,nThermal÷5;NSteps = nopt,NConfs,dt)

    flush(stdout)
    flush(stderr)

    CT = SW.ContinuousTimeMethod(tau,w_avg_estimate = E0_est,Hxx = CTSR.Hxx)
    initializer = SW.VariableTimePropagation(LinRange(10,tau,pre_init_steps),CT)

    @time results = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,estimate_w_avg=true,equilibration_steps=nThermal,initializer=initializer) for _ in 1:NRuns])

    results_en = [r.Energy[end] for r in results]

    results_sf = stack([r.StructureFactor[:,:,end] for r in results])

    flush(stdout)
    flush(stderr)

    return results_en,results_sf
end
##

sector_nums = [1,2,3,4,5,6,8,10]
muRange = LinRange(-0.1,0.99,30)

L = 20
Symms = [
    SW.TranslationalSymmetry([1,1],[1,-1]), #1
    SW.TranslationalSymmetry([2,2],[2,-2]), #2
    SW.TranslationalSymmetry([2,0],[0,4]), #3
    SW.TranslationalSymmetry([4,0],[0,4]), #4
    SW.TranslationalSymmetry([4,0],[0,4]), #5
    SW.TranslationalSymmetry([2,2],[0,4]), #6
    SW.TranslationalSymmetry([2,2],[2,-2]), #8
    SW.TranslationalSymmetry([4,0],[0,4]), #10
]
output_dir = "../../Data/EnergyScaling_S1_L$L"

function get_keys(file)
    h5open(file,"r") do f
        return keys(f)
    end
end

mkpath(output_dir)

println("start main loop")
flush(stdout)
flush(stderr)

for idx in eachindex(sector_nums,Symms)
    sector = sector_nums[idx]
    println("Sector $sector")
    flush(stdout)

    sector_dir = joinpath(output_dir, "sector_$sector")
    mkpath(sector_dir)

    Symm = Symms[idx]
    S = SW.get4x4PeriodicSpinConf(L,sector)

    for mu in muRange
        mu_file = joinpath(sector_dir, "mu_$(round(mu, digits=3)).h5")
        if isfile(mu_file)
            println("mu = $mu skipped")
            continue
        end
        try
            h5write(mu_file, "sector", sector)
            h5write(mu_file, "mu", mu)
            
        catch e
            @warn e
            continue
        end
        println("mu = $mu started")
        flush(stdout)

        SW.Random.seed!(1232)
        
        Nwalkers = 28*5
        nThermal = 500
        NSteps = 3000
        NRuns = 12
        nopt=1500
        NConfs=50
        pre_init_steps = 500

        if sector == 1
            Nwalkers *=20
            nThermal = 4000
            pre_init_steps = 1200
        end
        
        results_en, results_sf = getEnergy(S, mu, Symm; Nwalkers, NSteps, nThermal, NRuns, dt = 8e-3, NConfs, nopt,pre_init_steps)

        h5write(mu_file, "energies", results_en)
        h5write(mu_file, "structure_factors", results_sf)
        println("mu = $mu done")
    end
    println("Sector $sector done")
end

