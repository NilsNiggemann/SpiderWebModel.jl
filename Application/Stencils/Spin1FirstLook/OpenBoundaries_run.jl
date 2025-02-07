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

L = 24
SRoutfile = "../../Data/SR_open_L$L.h5"
mkpath(dirname(SRoutfile))

S = SW.stencilConfig(zeros(L,L),1,boundaryCondition = :open_soft)
ψG = SW.SimpleJastrowFunction(S)
SW.rand!(SW.get_params(ψG)) .*= 1e-6
ψG.v_ij .= SW.Symmetric(ψG.v_ij)

CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(0.0))
CTSR = SW.ContinuousTimeMethod(20*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

if !isfile(SRoutfile)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,50,ψG,1000,3e-3,SW.IterativeSRSolver();Nwalkers = 28*8,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=1_000,report_steps=10,outfile = SRoutfile)
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
outdir = "../../Data/open_L$(L)_3/"
rm(outdir,recursive=true,force=true)
mkdir(outdir)
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = w_avg_estimate,Hxx = CTSR.Hxx)
##
ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,28*100,4000,150,ψG,estimate_w_avg=true,equilibration_steps=2500,pre_equilibration_steps=10_000,outfile = "$outdir/$i.h5") for i in 1:12])
