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
S = SW.stencilConfig(zeros(24,24),1,boundaryCondition= :periodic) .= 2SW.periodicState6x6Condensate(24)
ψG = SW.SimpleJastrowFunction(S)
Symm = SW.TranslationalSymmetry([1,1],[1,-1])
ψGSymm = SW.symmetrize(ψG,Symm,S)
SW.rand!(ψGSymm,1e-3)
##
CTSR = SW.ContinuousTimeMethod(20*0.1,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.0))
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,20,ψGSymm,100,6e-3,SW.IterativeSRSolver();Nwalkers = 28*4,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=40_000,report_steps=5)
plotVarEn(stochReconfRes)
##
ObsRuns = [SW.measure_Sq_GFMC(S,CT,4000,3000,20,ψG,equilibration_steps=800,pre_equilibration_steps=10000) for _ in 1:1]

##
Sqs = [O.StructureFactor for O in ObsRuns]
with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    ax2 = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    SqMat = mean(Sqs)[:,:,20]
    SqErr = std(Sqs)[:,:,20]
    fittingCoefs = optimizeCoeffs(SqMat)
    
    Sq = SW.getSqCont(SqMat,cutoffEnd=0)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    heatmap!(ax,qx,qy, Sq.(qs)
    )
    SQFT(x) = SqFieldTheory(SVector(x),fittingCoefs)
    heatmap!(ax2,qx,qy,SQFT.(qs))
    fig

end