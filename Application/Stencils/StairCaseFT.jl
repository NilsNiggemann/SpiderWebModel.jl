import Pkg
Pkg.activate(joinpath(@__DIR__ ,"../"))
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
include("plottingUtils.jl")
##
S = SW.stencilConfig(zeros(16,16),1,boundaryCondition = :periodic) 
# S .= 4SW.getStairCase(size(S,1))
S .= SW.get4x4PeriodicState(16,6)
SW.plotApplPlaquettes(S)

##

SW.Random.seed!(1234)
mu = 0.8

CTSR = SW.ContinuousTimeMethod(10,Hxx = SW.Hxx_RK(mu))


nThermal = 1000

ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([2,2],[2,-2]),S)
SW.rand!(ψGSymm,1e-4)
w_avg_estimate = 0.1*length(S)

SRSteps = 800
##
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,30 ,ψGSymm,SRSteps,1e-2,SW.IterativeSRSolver();Nwalkers = 3*20,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 20,
reset = false,
# outfile = "tempSR/SR2.h5"
)
SW.get_params(ψG) .= stochReconfRes.params

w_avg_estimate = -stochReconfRes.E0[end]
plotVarEn(stochReconfRes,movavg = 20)
##
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = w_avg_estimate,Hxx = CTSR.Hxx)
##
results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20*20,5000,ψG;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:20])
##
plotEnergies(results,CT,normalize = false)

##
SqsGFMC = SW.getSqsGFMC(results,1:20:200)
##

with_theme(theme_SimpleTicks()) do 

    SqMat = dropmean(SqsGFMC,dims=4)[:,:,end,:]
    SqErr = dropstd(SqsGFMC,dims=4)[:,:,end,:]
    fittingCoefs = optimizeCoeffsAsym(SqMat)
    μ = mu
    fig = Figure(size = 120 .* (4,4),fontsize = 22)

    xticks = yticks = PiTicks([0,pi])
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks)

    ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = false,yticklabelsvisible = false)

    # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
    Sq = SW.getSqCont(SqMat)
    Sqerr = SW.getSqCont(SqErr)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    Sq_q = collect(Iterators.product(qx,qy))
    Sq_q = Sq.(Iterators.product(qx,qy))
    heatmap!(ax,qx,qy,Sq_q)
    
    SqFT = [AsymFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]

    heatmap!(axFT,qx,qy,SqFT)
    q_path(r,phi) = (r*cos(phi),r*sin(phi))
    qr = LinRange(0,.35pi,100)
    
    colors = (:red,:blue,:magenta)
    

    colorFT = :black
    colorGFMC = :red

    kpath = ["Γ","X","X'","Γ"]
    pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],500)
    kpointlabels = Makie.latexstring.(kpath)
    tRange = eachindex(p1)
    xygrid = [(x,y) for x in qx, y in qy]

    
    axPath = Axis(fig[2,1:2],ylabel = L"\mathcal{S}(\mathbf{q})" ,xlabel = L"\mathbf{q}" , xticks = (tRange[pointlabels],kpointlabels,),
    )
    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    

    p1_points = xygrid[p1_discrete]

    Sqcut = [Sq(x,y) for (x,y) in p1_points]
    Sqerrcut = [Sqerr(x,y) for (x,y) in p1_points]
    SqFT = [AsymFieldTheory(q,fittingCoefs) for q in p1_points]

    # SqFT = [AsymFieldTheory(q,1,10) for q in qpoints]
    scatter!(ax,p1_points,marker = '∘' ,color = colorGFMC,markersize = 15)
    scatterlines!(axFT,p1_points,color = colorFT,linestyle = :dash,marker = '●',markersize = 2)
    # tRange = SW.norm.(p1).^2
    scatterlines!(axPath,tRange,SqFT,color = colorFT,linestyle = :dash,marker = '●',markersize = 8)
    
    text!(axFT,Point(0,0),text="Γ",color = :white,align = (:center,:center))
    text!(axFT,Point(pi,0),text="X",color = :white,align = (:center,:center))
    text!(axFT,Point(0,pi),text="X'",color = :white,align = (:center,:center))

    scatter!(axPath,tRange,Sqcut,
    marker = '∘',markersize = 18,color = colorGFMC)
    errorbars!(axPath,tRange,Sqcut,Sqerrcut,color = colorGFMC,whiskerwidth = 6,linewidth=0.5)

    rowsize!(fig.layout,1,Relative(0.5))
    ylims!(axPath,0,6)
    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    # Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end
