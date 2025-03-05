using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

using ThreadPinning
ThreadPinning.pinthreads(:cores)

import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
##
# S = SW.stencilConfig(zeros(12,12),1,boundaryCondition= :periodic) .= 2SW.periodicState6x6Condensate(12)
S = SW.get4x4PeriodicSpinConf(20,1)
ψG = SW.SimpleJastrowFunction(S)
Symm = SW.TranslationalSymmetry([1,1],[1,-1])
ψGSymm = SW.symmetrize(ψG,Symm,S)
SW.rand!(ψGSymm,1e-5)
# SW.get_params(ψG) .= h5read("../../Data/GWF6x6/6x6_GF_24_mu0.h5","params")
##
CTSR = SW.ContinuousTimeMethod(20*0.1,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.2))
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,60,ψGSymm,2000,2e-2,SW.IterativeSRSolver();Nwalkers = 28*2,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=20_000,report_steps=25)
plotVarEn(stochReconfRes)
##
SW.get_params(ψG) .= stochReconfRes.params
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.18*length(S),Hxx = CTSR.Hxx)
# CT = SW.ContinuousTimeMethod(0.2,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.5))
# initializer = SW.CombinedInitializer(
#     SW.UnguidedWalkInitializer(10_000,0.9), 
#     SW.VariableTimePropagation(LinRange(8,0.1,500),CT)
# )
initializer = SW.CombinedInitializer(
    SW.UnguidedWalkInitializer(50_000,0.5), 
    # SW.StochasticResettingInitializer(LinRange(0.1,0.1,0),CT,100000.,S)
    # SW.StochasticResettingInitializer(exp10.(LinRange(0.9,log10(CT.τ),300)),CT,Inf,S)
    SW.StochasticResettingInitializer(LinRange(1,CT.τ,300),CT,40.,S)
)
# initializer =  SW.UnguidedWalkInitializer(50_000,0.9)
# initializer =  SW.StochasticResettingInitializer(exp10.(LinRange(0.8,log10(CT.τ),100)),CT,Inf,S)

steps_conv(Nsteps,CT) = round(Int,0.1/CT.τ*Nsteps)

DT = SW.DiscreteTimeMethod(0,21,4CT.w_avg_estimate)

# ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,100,steps_conv(1000,CT),ceil(Int,15/CT.τ),ψG;equilibration_steps=steps_conv(100,CT),initializer,estimate_w_avg=false) for _ in 1:50])
res = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,100,9000,ψG) for _ in 1:6])
equilib_plots(res,scatter_fraction=0.9)


# Sqs = [O.StructureFactor for O in ObsRuns[findall(x->minimum(x)>-90,energies)]]
# energies = [O.Energy for O in ObsRuns[findall(x->minimum(x)>-90,energies)]]

##
Sqs = [O.StructureFactor for O in ObsRuns]
energies = [O.Energy for O in ObsRuns]
filtered_indices = filter_outliers(energies)
println(length(filtered_indices))
# Sqs = setdiff(Sqs, Sqs[filtered_indices])
# energies = setdiff(energies, energies[filtered_indices])

Sqs = Sqs[filtered_indices]
energies = energies[filtered_indices]

with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 22,size = 500 .*(1.3,1.1))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    ax2 = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    SqMat = mean(Sqs)[:,:,1]
    SqErr = std(Sqs)[:,:,1]
    fittingCoefs = optimizeCoeffs(SqMat)
    
    Sq = SW.getSqCont(SqMat,cutoffEnd=0)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    heatmap!(ax,qx,qy, Sq)
    scatter!(ax,Point(pi/2,pi/2),marker = '∘',color = :black)
    SQFT(x...) = SqFieldTheory(SVector(x),fittingCoefs)
    heatmap!(ax2,qx,qy,SQFT)

    ax3 = Axis(fig[2, 1], xlabel = L"τ", ylabel = L"\mathcal{S}(q_i)",xticks=SimpleTicks(),yticks=SimpleTicks())


    tau_ax = CT.τ .* (eachindex(energies[1]) .-1)
    inds = [(1,2),(2,1),(5,4),[3,10]]

    for I in inds
        mean_Sq = mean(Sqs)[I...,:]
        err_Sq = std(Sqs)[I...,:]
        errlines!(ax3, tau_ax,mean_Sq, err_Sq)
    end
    mean_energy = mean(energies)
    err_energy = std(energies)

    ax4 = Axis(fig[2, 2], xlabel = L"τ", ylabel = L"E_0",xticks=SimpleTicks(),yticks=SimpleTicks())
    
    errlines!(ax4,tau_ax, mean_energy, err_energy)
    for (i, en) in enumerate(energies)
        lines!(ax4, tau_ax, en, linestyle = :dash, color = :black, linewidth = 0.4)

        textp = Point(rand(tau_ax), en[end])
        text!(ax4, textp, text = "$i", align = (:left, :center), color = :black,fontsize = 8,strokewidth=3,strokecolor = (:white,0.5))
        text!(ax4, textp, text = "$i", align = (:left, :center), color = :black,fontsize = 8)
    end
    lines!(ax4, [0, 0], [mean_energy[1] - err_energy[1], mean_energy[1] + err_energy[1]], color=:black)
    err_anot = string(round(err_energy[1],digits=5))
    text!(ax4, (0.5,0.8), text = L"%$(err_anot)", align = (:left, :bottom),space=:relative)
    ax5 = Axis(fig[3, :], xlabel = L"τ", ylabel = L"\mathcal{S}(0,2\pi/L)", xticks = SimpleTicks(), yticks = SimpleTicks())

    I = 1,2
    for (i, sq) in enumerate(Sqs)
        lines!(ax5, tau_ax, sq[I..., :], linestyle = :dash, color = :black, linewidth = 0.4)
        textp = Point(rand(tau_ax), sq[I..., end])
        text!(ax5, textp, text = "$i", align = (:left, :center), color = :black, fontsize = 8, strokewidth = 3, strokecolor = (:white, 0.5))
        text!(ax5, textp, text = "$i", align = (:left, :center), color = :black, fontsize = 8)
    end

    # xlims!(ax3,0,15)
    # xlims!(ax4,0,15)
    # xlims!(ax5,0,15)
    fig
end