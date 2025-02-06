import Pkg

import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(24,24),1,boundaryCondition= :periodic) .= 2SW.periodicState6x6Condensate(24)
ψG = SW.SimpleJastrowFunction(S)
Symm = SW.TranslationalSymmetry([6,0],[0,6])
ψGSymm = SW.symmetrize(ψG,Symm,S)
SW.rand!(ψGSymm,1e-5)
##
CTSR = SW.ContinuousTimeMethod(20*0.1,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.5))
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,80,ψGSymm,100,2e-3,SW.IterativeSRSolver();Nwalkers = 28*4,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=40_000,report_steps=5)
plotVarEn(stochReconfRes)