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
Confs = [SW.stencilConfig(float(S),1) for S in eachslice(h5read("../../Data/Spin1Confs/Confs_40_1.h5","Confs/1"),dims=3)]
Confs2 = [SW.stencilConfig(float(S),1) for S in eachslice(h5read("../../Data/Spin1Confs/Confs_40_2.h5","Confs/1"),dims=3)]
sConfs = vcat(Confs,Confs2)


sort!(sConfs,by=sum∘SW.getNPlaq,rev=true)
S = sConfs[1]
SW.plotApplPlaquettes(S)
##
ψG = SW.RKFunction()
# ψG = SW.SimpleJastrowFunction(S)
# ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([1,1],[1,-1]),S)
CTSR = SW.ContinuousTimeMethod(1.5,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(1))
# stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,200,ψGSymm,500,3e-4,SW.IterativeSRSolver();Nwalkers = 28*1,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=40_000,report_steps=1)
# plotVarEn(stochReconfRes)
##
CT = SW.ContinuousTimeMethod(2,Hxx = SW.Hxx_RK(1),w_avg_estimate=0.)

ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,28*2,2000,2,ψG,pre_equilibration_steps=500_000,equilibration_steps=100,estimate_w_avg=false,scatter_fraction=0.6) for i in 1:1*4])

##
with_theme(theme_PiTicks()) do 
    Sq = mean([OR.StructureFactor[:,:,1] for OR in ObsRuns])
    Sqerr = std([OR.StructureFactor[:,:,1] for OR in ObsRuns])
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig,ax,hm = heatmap(kx,ky,Sq,colormap = :viridis,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = (360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"stderr $$")
    heatmap!(ax2,kx,ky,Sqerr,colormap = :viridis,colorrange = extrema(Sq))
    Colorbar(fig[1:2,2],hm,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle")
    fig
end
##

