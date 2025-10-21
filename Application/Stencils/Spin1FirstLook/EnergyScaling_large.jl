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



function getEnergy(S,mu,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,pre_init_steps=800)
    ψG = SW.RKFunction()
    tau = 0.1 + 0.1*mu

    CTSR = SW.ContinuousTimeMethod(30*tau,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(mu))

    flush(stdout)
    flush(stderr)

    CT = SW.ContinuousTimeMethod(tau,w_avg_estimate = -0.05*length(S),Hxx = CTSR.Hxx)
    initializer = SW.UnguidedWalkInitializer(pre_init_steps,0.9)

    @time results = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,100,ψG,estimate_w_avg=true,equilibration_steps=nThermal,initializer=initializer) for _ in 1:NRuns])

    results_en = [r.Energy[end] for r in results]

    results_sf = stack([r.StructureFactor[:,:,end] for r in results])

    flush(stdout)
    flush(stderr)

    return results_en,results_sf
end
##

sector_nums = [1,2,3,4,5,6,8,10]
muRange = (0.9,)

L = 32
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
output_dir = ENV["MYSCRATCH"]*"/Spiderweb/4x4Comp/EnergyScaling_S1_L$L"

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
        
        Nwalkers = 28*100
        nThermal = 800
        NSteps = 8000
        NRuns = 12
        pre_init_steps = 100_000_000

        results_en, results_sf = getEnergy(S, mu, Symm; Nwalkers, NSteps, nThermal, NRuns, pre_init_steps)

        h5write(mu_file, "energies", results_en)
        h5write(mu_file, "structure_factors", results_sf)
        println("mu = $mu done")
    end
    println("Sector $sector done")
end

