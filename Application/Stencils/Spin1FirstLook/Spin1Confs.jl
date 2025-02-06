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
println("Threads: ",Threads.nthreads())
##
Confs = SW.constructGroundstatesSpin1(36,28*6,0.25,STotZero=true,TimeLimit = 200)

##

sConfs = SW.stencilConfig.(float.(Confs[1]),1)

sort!(sConfs,by=sum∘SW.getNPlaq,rev=true)
S = sConfs[1]
SW.plotApplPlaquettes(S)

##
ψG = SW.RKFunction()
CT = SW.ContinuousTimeMethod(4.0,Hxx = SW.Hxx_RK(1.0),w_avg_estimate=0.)

ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,28*10,3000,2,ψG,pre_equilibration_steps=10_000,equilibration_steps=1000,estimate_w_avg=false,scatter_fraction=0.95) for _ in 1:28])

##
with_theme(theme_PiTicks()) do 
    Sq = mean([OR.StructureFactor[:,:,1] for OR in ObsRuns])
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig,ax,hm = heatmap(kx,ky,Sq,colormap = :viridis,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = (360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,kx,ky,Sq,colormap = :viridis,colorrange = extrema(Sq))
    Colorbar(fig[1:2,2],hm,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle")
    fig
end
