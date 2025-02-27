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

L = 28
mu = 0.4
SRoutfile = "../../Data/SR_S0/SR_L_$(L)_mu_$(mu).h5"
mkpath(dirname(SRoutfile))

S = SW.stencilConfig(zeros(L,L),1,boundaryCondition = :periodic)
ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([1,1],[1,-1]),S)
SW.rand!(ψGSymm,1e-6)

CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(mu))
CTSR = SW.ContinuousTimeMethod(30*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

if !isfile(SRoutfile)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,60,ψGSymm,4000,8e-3,SW.IterativeSRSolver();Nwalkers = 28*3,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=1_000,report_steps=25,outfile = SRoutfile)
    # h5write(SRoutfile,"E0",stochReconfRes.E0)
    # h5write(SRoutfile,"ΔE",stochReconfRes.ΔE)
    # h5write(SRoutfile,"params",stochReconfRes.params)
end
#___________Open Boundaries_______________________

##
idxzero,w_avg_estimate = h5open(SRoutfile,"r") do f
    E0 = read(f["E0"])
    idx0 = findfirst(iszero,E0)
    if isnothing(idx0)
        return lastindex(E0),E0[end]
    end
    return idx0,E0[idx0]
end

SW.get_params(ψG) .= h5read(SRoutfile,"params_steps")[:,idxzero]
outdir = "../../Data/S0periodic_L$(L)/mu=$mu"
rm(outdir,recursive=true,force=true)
mkpath(outdir)
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = w_avg_estimate,Hxx = CTSR.Hxx)
##
initializer = SW.CombinedInitializer(
    SW.UnguidedWalkInitializer(50_000,0.4), 
    SW.StochasticResettingInitializer(exp10.(LinRange(1,log10(CT.τ),1000)),CT,200.,S)
)
ObsRuns = [SW.measure_Sq_GFMC(S,CT,28*300,4000,80,ψG;estimate_w_avg=true,equilibration_steps=2000,outfile = "$outdir/$i.h5",initializer) for i in 1:12]
# ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,28*10,2000,80,ψG;estimate_w_avg=true,equilibration_steps=1500,pre_equilibration_steps=10_000,outfile = "$outdir/$i.h5",initializer) for i in 1:4])


