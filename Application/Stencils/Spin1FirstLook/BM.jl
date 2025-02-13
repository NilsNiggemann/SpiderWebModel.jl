using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

using ThreadPinning
ThreadPinning.pinthreads(:cores)

using Pkg
Pkg.activate(joinpath(@__DIR__,"../../.."))
import SpiderWebModel as SW
using ThreadPinning
pinthreads(:cores)

S = SW.stencilConfig(zeros(30,30),1,boundaryCondition=:periodic)
S .= 2*SW.periodicState6x6Condensate(size(S,1))
# S = SW.stencilConfig(zeros(20,20),1,boundaryCondition=:periodic)
ψG = SW.SimpleJastrowFunction(S)
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.1))
a = SW.startManyWalkerGFMC(S,CT,10*14,10,ψG);
# ψGSymm = SW.getNonSymmetric(ψG)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([6,0],[0,6]),S)
SW.rand!(ψGSymm,1e-5)
##
CTSR = SW.ContinuousTimeMethod(5,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.1))
##
SW.Random.seed!(1234)
@profview @time SW.stochastic_reconfiguration(S,CTSR,50,ψGSymm,10,1e-4,SW.IterativeSRSolver();Nwalkers = 28*4,equilibration_steps=0,report_steps=5)

# @profview @time a = SW.startManyWalkerGFMC(S,CT,4000*14,10,ψG);
# @time a = SW.startManyWalkerGFMC(S,CT,50*14,20,ψG);

nothing