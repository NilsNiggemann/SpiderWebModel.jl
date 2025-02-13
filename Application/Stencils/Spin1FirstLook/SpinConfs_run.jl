using LinearAlgebra
BLAS.set_num_threads(Threads.nthreads())
using ThreadPinning
ThreadPinning.pinthreads(:cores)
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
println("Threads: ",Threads.nthreads())
##
Confs = [SW.stencilConfig(float(S),1) for S in eachslice(h5read("../../Data/Spin1Confs/Confs_40_1.h5","Confs/15"),dims=3)]
Confs2 = [SW.stencilConfig(float(S),1) for S in eachslice(h5read("../../Data/Spin1Confs/Confs_40_2.h5","Confs/15"),dims=3)]
sConfs = vcat(Confs,Confs2)


sort!(sConfs,by=sum∘SW.getNPlaq,rev=true)
S = sConfs[1]
# ψG = SW.RKFunction()
ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([1,1],[1,-1]),S)
CTSR = SW.ContinuousTimeMethod(3,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.6))
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,400,ψGSymm,2000,2e-4,SW.IterativeSRSolver();Nwalkers = 28*3,rel_tolerance=1e-8,equilibration_steps=500,pre_equilibration_steps=40_000,report_steps=1,outfile = "../../Data/Spin1Confs/SR_40_1.h5")