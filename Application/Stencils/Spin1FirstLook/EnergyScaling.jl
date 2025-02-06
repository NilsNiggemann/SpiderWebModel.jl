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
        idx > 5 && error("Optimization failed")
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

function getEnergies(S,mus,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,dt = 5e-3,NConfs = 20,nopt1 = 100,nopt2 = 30)
    ψG = initializeGWF(S,Symmetry)
    ens = zeros(length(mus),NRuns)
    for (i,mu) in enumerate(mus)
        println("mu = $mu")

        tau = 0.1 + 0.1*mu
        CTSR = SW.ContinuousTimeMethod(50*tau,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(mu))

        nopt = i == 1 ? nopt1 : nopt2
        @time E0_est = optimizeWF!(ψG,S,CTSR,nThermal÷5;NSteps = nopt,NConfs,dt)
        flush(stdout)
        flush(stderr)

        CT = SW.ContinuousTimeMethod(tau,w_avg_estimate = E0_est,Hxx = CTSR.Hxx)

        @time results_en = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,estimate_w_avg=true,equilibration_steps=nThermal,pre_equilibration_steps=100*nThermal).Energy for _ in 1:NRuns])
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

println("start main loop")
flush(stdout)
flush(stderr)

# for (sector,Symm) in zip(sector_nums,Symms)
for idx in eachindex(sector_nums,Symms)
    sector = sector_nums[idx]
    print("Sector $sector")

    present_keys = get_keys(outfile_EnergyScaling)
    if string(sector) in present_keys
        println(" skipped")
        flush(stdout)
        continue
    end
    println(" started")
    flush(stdout)
    
    Symm = Symms[idx]
    S = SW.get4x4PeriodicSpinConf(L,sector)
    h5write(outfile_EnergyScaling,"$sector/conf", Array(S))

    SW.Random.seed!(1232)
    
    Nwalkers = 28*20
    nThermal = 2000

    if sector == 1
        Nwalkers *=12
        nThermal = 4000
    end
    e_i = getEnergies(S,muRange,Symm;Nwalkers,NSteps = 3000,nThermal,NRuns=12,dt = 3e-3,NConfs = 50,nopt1 = 700,nopt2 = 600)
    
    h5write(outfile_EnergyScaling,"$sector/energy",e_i)
    println("Sector $sector done")
end
exit()
##
using SpiderWebModel.HDF5
using MakieHelpers
using CairoMakie
function plot_energies(outfile)
    fig = Figure(size = (700, 500),fontsize = 22)
    

    toprow = fig[1,1:4] = GridLayout()
    energyFig = fig[2:4,1:4] = GridLayout()
    # Top row for configurations

    h5open(outfile, "r") do f
        
        sectors = sort(filter(!=("muRange"),keys(f)),by=x->parse(Int,x))
        
        colors = Makie.ColorSchemes.distinguishable_colors(length(sectors), lchoices = LinRange(0,80,100))

        for (i, sector) in enumerate(sectors)
            conf = read(f["$sector/conf"])
            if sector == "6"
                conf .= 4SW.getStairCase(size(conf,1)) #plot the staircase state, which is in the same sector
            end
            S = SW.stencilConfig(conf[1:8,1:8],1)
            toprow[1,i] = ax_conf = Axis(fig, title=L"%$i";SW.getConfigAxis(S)...,
            yticklabelsvisible=false,
            xticklabelsvisible=false,
            spinecolors(colors[i])...,
            spinewidth = 4
            )
            SW.plotSpinConfig!(ax_conf, S)
            colgap!(toprow, 0)
        end

        # Bottom row for energies
        energyFig[1,1:8] = ax = with_theme(theme_SimpleTicks()) do
            Axis(fig, xlabel=L"\mu", 
            ylabel=L"E_0/N_{\text{sites}}"
            # ylabel=L"E_0/(N_{\text{sites}}(1-\mu))"
            )
        end
        L = 20
        muRange = f["muRange"][:]
        for (i, sector) in enumerate(sectors)
            "energy" in keys(f[sector]) || continue
            Nsites = length(f[sector*"/conf"])
            energies = read(f["$sector/energy"])
            mean_energies = dropmean(energies, dims=2) ./Nsites# ./ (1 .-muRange)
            std_energies = dropstd(energies, dims=2) ./Nsites# ./ (1 .-muRange)

            errlines!(ax, muRange, mean_energies, std_energies, label=L"%$i", markersize=0.1,color = colors[i],linewidth = 2)
        end
        axislegend(ax, position=:lt, merge=true)
    end
    fig
end
outfile_EnergyScaling = "../../Data/energy_mu_S1.h5"

plot_energies(outfile_EnergyScaling)
##

function recoverCorrupt(fname,newfile,keys)
    try
        muRange = h5read(fname,"muRange")
        h5write(newfile,"muRange",muRange)
        
    catch
    end
    for key in keys, key2 in ("conf","energy")
        try
            data = h5read(fname,"$key/$key2")
            h5write(newfile,"$key/$key2",data)
        catch

        end
    end

end
recoverCorrupt(outfile_EnergyScaling,"../../Data/energy_mu_S1_recovered.h5",["1","2","3","4","5","6","8","10"])
