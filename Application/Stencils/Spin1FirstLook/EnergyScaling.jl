using ThreadPinning
ThreadPinning.pinthreads(:cores)
import Pkg
cd(@__DIR__)
Pkg.activate("../../")
# Pkg.instantiate()
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using SpiderWebModel.HDF5
include("../plottingUtils.jl")
println("Threads: ",Threads.nthreads())
##
#___________Periodic Boundaries_______________________

function initializeGWF(S,Symmetry) 
    ψG = SW.SimpleJastrowFunction(S)
    ψGSymm = SW.symmetrize(ψG,Symmetry,S)
    SW.rand!(ψGSymm,1e-5)
    return ψGSymm
end

function _optimizeWF(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,NConfs,ψG,NSteps,dt,SW.IterativeSRSolver();Nwalkers = 20,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=25,kwargs...)
    if showplot
        display(plotVarEn(stochReconfRes))
    end
    return stochReconfRes
end

function optimizeWF!(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = _optimizeWF(ψG,S,CTSR,nThermal;dt = dt,NSteps = NSteps,NConfs = NConfs,showplot=showplot,kwargs...)

    E0 = stochReconfRes.E0
    idx = 0
    while E0[end] > E0[begin]
        idx += 1
        idx > 5 && error("Optimization failed")
        dt *= 0.5
        NSteps = NSteps *2
        stochReconfRes = _optimizeWF(ψG,S,CTSR,nThermal;dt = dt,NSteps = NSteps,NConfs = NConfs,showplot=showplot,kwargs...)
        E0 = stochReconfRes.E0
    end
    SW.get_params(ψG) .= stochReconfRes.params

    return stochReconfRes.E0[end]
end
##

function getEnergies(S,mus,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,dt = 5e-3,NConfs = 20,nopt1 = 100,nopt2 = 30)
    ψG = initializeGWF(S,Symmetry)
    ens = zeros(length(mus),NRuns)
    for (i,mu) in enumerate(mus)
        println("mu = $mu")

        tau = 0.1 + 0.1*mu
        CTSR = SW.ContinuousTimeMethod(30*tau,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(mu))

        nopt = i == 1 ? nopt1 : nopt2
        @time E0_est = optimizeWF!(ψG,S,CTSR,nThermal÷5;NSteps = nopt,NConfs,dt)
        flush(stdout)
        flush(stderr)

        CT = SW.ContinuousTimeMethod(tau,w_avg_estimate = E0_est,Hxx = CTSR.Hxx)

        @time results_en = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,estimate_w_avg=true,equilibration_steps=nThermal).Energy for _ in 1:NRuns])
        flush(stdout)
        flush(stderr)
        # @time results_en = SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,equilibration_steps=NSteps÷5, estimate_w_avg=true)
        # return results_en
        # display(plotEnergies(results,CT))
        en = stack(results_en)[end,:]
        ens[i,:] .= en

    end
    return ens
end
##

sector_nums = [1,2,3,4,5,6,8,10]
muRange = reverse(LinRange(-0.1,0.99,30))

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
outfile_EnergyScaling = "../../Data/energy_mu_S1.h5"

##

# outfile_EnergyScaling = "../../Data/energy_mu_S1.h5"

function get_keys(file)
    h5open(file,"r") do f
        return keys(f)
    end
end


h5open(outfile_EnergyScaling,"cw") do f
    
    pres_keys = keys(f)
    if "muRange" in pres_keys
        if f["muRange"][:] != Array(muRange)
            rm(outfile_EnergyScaling)
            f["muRange"] = Array(muRange)
        end
    else    
        f["muRange"] = Array(muRange)
    end
end

mkpath(dirname(outfile_EnergyScaling))


# for (sector,Symm) in zip(sector_nums,Symms)
for idx in eachindex(sector_nums,Symms)
    sector = sector_nums[idx]
    flush(stdout)
    present_keys = get_keys(outfile_EnergyScaling)
    string(sector) in present_keys && continue
    println("Sector $sector")

    Symm = Symms[idx]
    S = SW.get4x4PeriodicSpinConf(L,sector)
    h5write(outfile_EnergyScaling,"$sector/conf", Array(S))

    SW.Random.seed!(1232)
    
    e_i = getEnergies(S,muRange,Symm,Nwalkers = 1000,NSteps = 3000,nThermal = 2000,NRuns=30,dt = 4e-4,NConfs = 20,nopt1 = 400,nopt2 = 700)
    
    h5write(outfile_EnergyScaling,"$sector/energy",e_i)
    println("Sector $sector done")
end
exit()
##
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
function plot_energies(outfile)
                
    fig = Figure()
    ax = with_theme(theme_SimpleTicks()) do
        Axis(fig[1, 1], xlabel=L"\mu", ylabel=L"E_0/(N_{\text{sites}}(1-\mu))")
    end
    L = 20
    h5open(outfile, "r") do f
        muRange = f["muRange"][:]
        for sector in keys(f)
            sector == "muRange" && continue
            "energy" in keys(f[sector]) || continue
            energies = read(f["$sector/energy"])
            mean_energies = dropmean(energies, dims=2) ./ L^2 ./ (1 .-muRange)
            std_energies = dropstd(energies, dims=2) ./ L^2 ./ (1 .-muRange)

            # errorbars!(ax, muRange, mean_energies, yerr=std_energies)
            errlines!(ax, muRange, mean_energies, std_energies,label = L"%$sector",markersize = 5)
        end
    end
    axislegend(ax, position = :lt)
    fig
end
# outfile_EnergyScaling = "../../Data/energy_mu_S1_2.h5"


plot_energies(outfile_EnergyScaling)