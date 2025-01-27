import Pkg
cd(@__DIR__)
Pkg.activate("../../")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using SpiderWebModel.HDF5
include("../plottingUtils.jl")
##

L = 20
SRoutfile = "../../Data/SR_open_L$L.h5"
mkpath(dirname(SRoutfile))

S = SW.stencilConfig(zeros(20,20),1,boundaryCondition = :open_soft)
ψG = SW.SimpleJastrowFunction(S)
SW.rand!(SW.get_params(ψG)) .*= 1e-6
ψG.v_ij .= SW.Symmetric(ψG.v_ij)

CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(0.0))
CTSR = SW.ContinuousTimeMethod(50*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

if !isfile(SRoutfile)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,150,ψG,1000,3e-3,SW.IterativeSRSolver();Nwalkers = 40,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=1_000,report_steps=10,outfile = SRoutfile)
    # h5write(SRoutfile,"E0",stochReconfRes.E0)
    # h5write(SRoutfile,"ΔE",stochReconfRes.ΔE)
    # h5write(SRoutfile,"params",stochReconfRes.params)
end
#___________Open Boundaries_______________________

##
SW.get_params(ψG) .= h5read(SRoutfile,"params_steps")[:,end]
rm("../../Data/open_L$L/",recursive=true,force=true)
mkdir("../../Data/open_L$L/")
ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,500,3000,150,ψG,estimate_w_avg=true,equilibration_steps=1500,pre_equilibration_steps=10_000,outfile = "../../Data/open_L$L/$i.h5") for i in 1:20])
