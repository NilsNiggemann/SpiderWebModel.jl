cd(@__DIR__)
using Pkg
Pkg.activate("../..")
using LinearAlgebra
println("Threads: ",Threads.nthreads())
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

using ThreadPinning
ThreadPinning.pinthreads(:cores)

using CairoMakie, MakieHelpers,Statistics
import SpiderWebModel as SW

include("../plottingUtils.jl")
SW.initGurobi()
##

# SW.Random.seed!(1233)
SW.Random.seed!(100)
L = 40
@time sols = SW.constructGroundstatesSpin1(L, 500, 0.2,boundaryCondition = :periodic,STot=(0,0),S_staggered = 0,MI=zeros(L),M_ = zeros(L),M_diag = zeros(L),M_anti = zeros(L),progress = true,TimeLimit = 20)
##
confs = [SW.stencilConfig(float(x),1,boundaryCondition = :periodic) for x in sols.solutions]
##
sort!(confs,by = SW.NPlaquettes,rev=true)
S = confs[1]
SW.plotApplPlaquettes(S)
##
function findMaxFlipConf(S;Nwalkers=200,NSteps =1000,kwargs...)
    CT = SW.ContinuousTimeMethod(0.5,w_avg_estimate = 0.,Hxx = SW.Hxx_RK(-10))
    ψG = SW.PlaquetteNumberGuidingFunction(1.)
    res = SW.startManyWalkerGFMC(S,CT,Nwalkers,NSteps,ψG,pre_equilibration_steps=800_000,scatter_fraction=0.9;kwargs...)

    AllConfs = [SW.stencilConfig(Array(x),1) for x in eachslice(res.SaveConfigs,dims = (3,4))]

    maxConf = argmax(SW.NPlaquettes,AllConfs)
    return maxConf

end
S = findMaxFlipConf(S;Nwalkers=500,NSteps =2000)
SW.plotApplPlaquettes(S)
##
maxConfs = [findMaxFlipConf(confs[i];Nwalkers=100,NSteps =200) for i in 1:10]
##
let 
    fig = Figure(resolution = 0.8 .*(1200, 900))
    n = 10
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)
    for i in 1:n
        row = div(i-1, ncols) + 1
        col = mod(i-1, ncols) + 1
        ax = Axis(fig[row, col]; SW.getConfigAxis(S)...,xticklabelsvisible = false,yticklabelsvisible = false)
        SW.plotApplPlaquettes!(ax, maxConfs[i],markersize = 4)
    end
    fig
end
##
# CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = w_avg_estimate,Hxx = CTSR.Hxx)
mu = 1
# CT = SW.ContinuousTimeMethod(0.5,w_avg_estimate = 0.,Hxx = SW.Hxx_RK(mu))
CT = SW.ContinuousTimeMethod(0.5,w_avg_estimate = 0.,Hxx = SW.Hxx_RK(mu))
ψG = SW.RKFunction()
nThermal = 1000
##
# initializer = SW.CombinedInitializer(
#     SW.UnguidedWalkInitializer(100_000_000, 0.9),
#     # SW.StochasticResettingInitializer(LinRange(2, CT.τ, 400), CT_init, Inf, confs[1])
#     SW.StochasticResettingInitializer(LinRange(1, CT.τ, 400), CT, Inf, S)
# )
ObsRuns = [SW.measure_Sq_GFMC(maxConfs[1], CT, 28, 2000, 2, ψG, estimate_w_avg=false,scatter_fraction=1.,pre_equilibration_steps=50_000_000)]

# results = [SW.startManyWalkerGFMC(confs[1],CT,8,20000,ψG,scatter_fraction=0.9,pre_equilibration_steps=800_000)]
SqsGFMC = [O.StructureFactor for O in ObsRuns]
# SqsGFMC = [SW.getSqsGFMC(results,1:2)[:,:,1]]
##

function plotSq(SqsGFMC)
    with_theme(theme_SimpleTicks()) do 
        SqMat = mean(SqsGFMC)[:,:,end]
        # SqMat[end÷2+1,end÷2+1] = NaN
        SqErr = std(SqsGFMC)[:,:,end]
        fittingCoefs = optimizeCoeffs(SqMat)
        μ = mu
        fig = Figure(size = 120 .* (4,4),fontsize = 22)

        xticks = yticks = PiTicks([0,pi])
        axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks)

        ax = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1;xticks, yticks,ylabelvisible = false,yticklabelsvisible = false)

        # ax2 = Axis(fig[2,1:2],xlabel = L"|\mathbf{q}|^2",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"μ= %$μ")
        Sq = SW.getSqCont(SqMat,cutoffEnd=0)
        Sqerr = SW.getSqCont(SqErr)
        qx = qy = trueMomenta(-0.5pi,1.5pi,size(SqMat,1))
        heatmap!(ax,qx,qy,Sq)
        
        SqFT = [SqFieldTheory(x,y,fittingCoefs...) for x in qx, y in qy]
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
        SqFT = [SqFieldTheory(q,fittingCoefs) for q in p1_points]

        # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
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
        # ylims!(axPath,0,6)

        fig
    end    
end
plotSq(SqsGFMC)
##
SW.Random.seed!(1234)
mu = 0.8

CTSR = SW.ContinuousTimeMethod(0.5,Hxx = SW.Hxx_RK(mu))


nThermal = 100

ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([1,1],[1,-1]),S)
SW.rand!(ψGSymm,1e-4)
w_avg_estimate = 0.1*length(S)
##
using MKL
@time stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,80 ,ψGSymm,1000,5e-2,SW.IterativeSRSolver();Nwalkers = 80*2,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
report_steps = 20,
reset = false,
# outfile = "tempSR/SR2.h5"
)

w_avg_estimate = -stochReconfRes.E0[end]
plotVarEn(stochReconfRes,movavg = 20)
##
SW.get_params(ψG) .= stochReconfRes.params

CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = w_avg_estimate,Hxx = CTSR.Hxx)

initializer = SW.CombinedInitializer(
    SW.UnguidedWalkInitializer(80_000_000, 0.9),
    SW.StochasticResettingInitializer(LinRange(1, CT.τ, 300), CT, Inf, S)
)

ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S, CT, 1000, 5000, 40, ψG, initializer=initializer,equilibration_steps = 1500,estimate_w_avg=true) for _ in 1:6])

SqsGFMC = [O.StructureFactor for O in ObsRuns]
En = [O.Energy for O in ObsRuns]
ids = filter_outliers(En)
En = En[ids]
SqsGFMC = SqsGFMC[ids]
errlines(mean(En),std(En))
##
plotSq(SqsGFMC)
##

results = [SW.startManyWalkerGFMC(S,CT,1000,5000,ψG;initializer)]
##

SW.plotApplPlaquettes(copy(S) .= results[1].SaveConfigs[:,:,end,end-40])

##
S0 = copy(S) .= 0
initializer = SW.CombinedInitializer(
    SW.UnguidedWalkInitializer(80_000_000, 1.),
    SW.StochasticResettingInitializer(LinRange(1, CT.τ, 200), CT, Inf, S0)
)
resultsS0 = [SW.startManyWalkerGFMC(S0,CT,1000,5000,ψG;initializer)]
##

unique_configs = Set(Array(x) for x in eachslice(results[1].SaveConfigs,dims = (3,4)))
unique_configs_S0 = Set(Array(x) for x in eachslice(resultsS0[1].SaveConfigs,dims = (3,4)))

inters = intersect(unique_configs,unique_configs_S0)

##
SW.plotApplPlaquettes(copy(S) .= resultsS0[1].SaveConfigs[:,:,end,end-1000])
##
ObsRunsS0 = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S0, CT, 100, 2000, 40, ψG, scatter_fraction=1.0, initializer=initializer,equilibration_steps = 1000,estimate_w_avg=true) for _ in 1:10])

SqsGFMC_S0 = [O.StructureFactor for O in ObsRunsS0]
En_S0 = [O.Energy for O in ObsRunsS0]
errlines(mean(En_S0),std(En_S0))
##
plotSq(SqsGFMC_S0)