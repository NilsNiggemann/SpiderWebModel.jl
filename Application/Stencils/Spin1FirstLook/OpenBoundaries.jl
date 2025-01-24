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

#___________Open Boundaries_______________________

S = SW.stencilConfig(zeros(20,20),1,boundaryCondition = :open_soft)

GFMCfiles = readdir("../../Data/open_L20/",join=true)

SqsGFMC_direct = stack([h5read(file,"StructureFactor") for file in GFMCfiles])

with_theme(theme_SimpleTicks()) do 
    SqMat = SW.expand_Sq(dropmean(SqsGFMC_direct,dims=4))[:,:,50]
    SqErr = SW.expand_Sq(dropstd(SqsGFMC_direct,dims=4))[:,:,50]
    fittingCoefs = optimizeCoeffs(SqMat)
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
    # yscale = log10
    )
    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    

    p1_points = xygrid[p1_discrete]

    Sqcut = [Sq(x,y) for (x,y) in p1_points]
    Sqerrcut = [Sqerr(x,y) for (x,y) in p1_points]
    SqFT = [SqFieldTheory(q,fittingCoefs) for q in p1_points]

    # SqFT = [AsymFieldTheory(q,1,10) for q in qpoints]
    scatter!(ax,p1_points,marker = '∘' ,color = colorGFMC,markersize = 15)
    scatterlines!(axFT,p1_points,color = colorFT,linestyle = :dash,marker = '●',markersize = 2)
    # tRange = SW.norm.(p1).^2
    scatterlines!(axPath,tRange,SqFT,color = colorFT,linestyle = :dash,marker = '●',markersize = 8)
    
    text!(axFT,Point(0,0),text="Γ",color = :white,align = (:center,:center))
    text!(axFT,Point(pi,0),text="X",color = :white,align = (:center,:center))
    text!(axFT,Point(0,pi),text="X'",color = :white,align = (:center,:center))

    scatter!(axPath,tRange,Sqcut.+1e-30,
    marker = '∘',markersize = 18,color = colorGFMC)
    errorbars!(axPath,tRange,Sqcut.+1e-30,Sqerrcut,color = colorGFMC,whiskerwidth = 6,linewidth=0.5)

    rowsize!(fig.layout,1,Relative(0.5))
    # ylims!(axPath,0.01,nothing)
    # ylims!(axPath,0,1.2)
    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    # Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end


##
function initializeGWF(S,Symmetry) 
    ψG = SW.SimpleJastrowFunction(S)
    ψGSymm = SW.symmetrize(ψG,Symmetry,S)
    SW.rand!(ψGSymm,1e-5)
    return ψGSymm
end

function optimizeWF!(ψG,S,CTSR,nThermal;dt = 1e-4,NSteps = 100,NConfs = 10,showplot=false,kwargs...)
    stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,NConfs,ψG,NSteps,dt,SW.IterativeSRSolver();Nwalkers = 20,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=25,kwargs...)
    SW.get_params(ψG) .= stochReconfRes.params
    if showplot
        display(plotVarEn(stochReconfRes))
    end
    return ψG
end

##

function getEnergies(S,mus,Symmetry;Nwalkers = 20 * 3,NSteps = 2000,nThermal = 1000,NRuns=20,dt = 5e-3,NConfs = 20,nopt1 = 100,nopt2 = 30)
    ψG = initializeGWF(S,Symmetry)
    ens = zeros(length(mus),NRuns)
    for (i,mu) in enumerate(mus)
        println("mu = $mu")
        CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(mu))
        CTSR = SW.ContinuousTimeMethod(10*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

        nopt = i == 1 ? nopt1 : nopt2
        @time optimizeWF!(ψG,S,CTSR,nThermal;NSteps = nopt,NConfs,dt)
        @time results_en = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,estimate_w_avg=false).Energy for _ in 1:NRuns])
        # @time results_en = SW.measure_Sq_GFMC(S,CT,Nwalkers,NSteps,50,ψG.psi,equilibration_steps=NSteps÷5, estimate_w_avg=true)
        # return results_en
        # display(plotEnergies(results,CT))
        en = stack(results_en)[end,:]
        ens[i,:] .= en

    end
    return ens
end

function mainRun(Ls,muRange,Symmetries)
    run(i) = getEnergies( getSectorConfig(Ls[i],i),muRange,Symmetries[i])

    allE = fetch.([Threads.@spawn run(i) for i in 1:3])

    return allE
end
##

muRange = reverse(LinRange(0,0.99,20))

Ls = [20,24,20,18]
Symms = [
    SW.TranslationalSymmetry([-2,2],[2,2]),
    SW.TranslationalSymmetry([6,0],[0,6]),
    SW.TranslationalSymmetry([5,1],[0,4]),
    SW.TranslationalSymmetry([3,-1],[0,6])
]
ens1 = getEnergies( getSectorConfig(Ls[1],1),muRange,Symms[1],NConfs=20,dt = 2e-4)
##
ens2 = getEnergies( getSectorConfig(Ls[2],2),muRange,Symms[2],NConfs=20,dt = 2e-4)
##
ens3 = getEnergies( getSectorConfig(Ls[3],3),muRange,Symms[3],NConfs=20,dt = 2e-4)
##
ens4 = getEnergies( getSectorConfig(Ls[4],4),muRange,Symms[4],NConfs=20,dt = 2e-4)
# ens2 = mainRun([20,20,20,20],muRange)
##
h5write("../../Data/energy_mu_S12.h5","muRange", Array(muRange))
h5write("../../Data/energy_mu_S12.h5","Sector1", Array(getSectorConfig(Ls[1],1)))
h5write("../../Data/energy_mu_S12.h5","energy1",ens1)
h5write("../../Data/energy_mu_S12.h5","Sector2", Array(getSectorConfig(Ls[2],2)))
h5write("../../Data/energy_mu_S12.h5","energy2",ens2)
h5write("../../Data/energy_mu_S12.h5","Sector3", Array(getSectorConfig(Ls[3],3)))
h5write("../../Data/energy_mu_S12.h5","energy3",ens3)
h5write("../../Data/energy_mu_S12.h5","Sector4", Array(getSectorConfig(Ls[4],4)))
h5write("../../Data/energy_mu_S12.h5","energy4",ens4)
##

S = getSectorConfig(32,1)
ψGSymm = initializeGWF(S,SW.TranslationalSymmetry([-2,2],[2,2]))
ψG = ψGSymm.psi
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.1*length(S),Hxx = SW.Hxx_RK(0.0))
CTSR = SW.ContinuousTimeMethod(30*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

# @time optimizeWF!(ψG,S,CTSR,1000;NSteps = 1000,NConfs=200,dt=1e-3,showplot=true,report_steps = 4)
# @time optimizeWF!(ψGSymm,S,CTSR,100;NSteps = 100,NConfs=200,dt=1e-3,showplot=true,report_steps = 4)

stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,300,ψGSymm,700,1e-3,SW.IterativeSRSolver();Nwalkers = 20,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=40_000,report_steps=10,outfile = "../../Data/Stoch_reconf_S12_Stair.h5")
SW.get_params(ψG) .= stochReconfRes.params
plotVarEn(stochReconfRes)
##
SW.Random.seed!(1232)
ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,300,5000,150,equilibration_steps=2000,ψG,estimate_w_avg=true,outfile = "../../Data/obsS12_staircase/obsS12_staircase_mu0_$i.h5") for i in 1:12])
##
# ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,20*10,2000,150,ψG,equilibration_steps=1000,estimate_w_avg=true) for i in 1:12])
results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,200,1000,ψG,equilibration_steps=2000) for _ in 1:20])
##
SW.Random.seed!(1232)
ObsRuns = fetch.([Threads.@spawn SW.measure_Sq_GFMC(S,CT,300,500,150,equilibration_steps=2000,ψG,estimate_w_avg=true,) for i in 1:12])
##
en_direct = stack([ObsRuns.Energy for ObsRuns in ObsRuns])

plotEnergies(results,CT,nThermal=1,τ=10,normalize=false)
# plotEnergies!(ObsRuns,CT,nThermal=1,τ=10,normalize=false,color = :red)
errlines!((0:size(en_direct,1)-1).*CT.τ,dropmean(en_direct,dims=2),dropstd(en_direct,dims=2))
current_figure()
##
Sqs = stack([SW.expand_Sq(O.StructureFactor) for O in ObsRuns])
with_theme(theme_PiTicks()) do
    fig = Figure(fontsize = 22,size = (400,320))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks((0,pi)),yticks = PiTicks((0,pi)))
    # ax2 = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",aspect=1)
    SqMat = dropmean(Sqs,dims=4)[:,:,end]
    SqError = dropstd(Sqs,dims=4)[:,:,end]

    # fittingCoefs = optimizeCoeffs(SqMat)
    Sq = SW.getSqCont(SqError)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    hm = heatmap!(ax,qx,qy, Sq.(qs))
    # hm = heatmap!(ax,qx,qy, Sq.(qs) .+1e-10,colorscale = log10,colorrange = (0.01,maximum(SqMat)))
    Colorbar(fig[1,2],hm,label = L"\langle S^z(\mathbf{q})^*  S^z(\mathbf{q})\rangle")
    # SQFT(x) = SqFieldTheory(SVector(x),fittingCoefs)
    # heatmap!(ax2,qx,qy,SQFT.(qs))
    fig
end

##
S = getSectorConfig(24,1)
ψG = SW.RKFunction()

CT = SW.ContinuousTimeMethod(0.5,w_avg_estimate = 0.,Hxx = SW.Hxx_RK(0.0))

resultsRK = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20*1,8000,ψG,pre_equilibration_steps = 100000) for _ in 1:20]) 
##
SqsRK = SW.getSqsGFMC(resultsRK,1:100)
with_theme(theme_PiTicks()) do
    fig = Figure(fontsize = 22,size = (400,320))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect=1,xticks = PiTicks((0,pi)),yticks = PiTicks((0,pi)))
    SqMat = dropmean(SqsRK,dims=4)[:,:,end]
    SqError = dropstd(SqsRK,dims=4)[:,:,end]

    Sq = SW.getSqCont(SqMat)
    qx = qy = trueMomenta(-0.5pi,1.5pi,size(S,1))
    qs = Iterators.product(qx,qy)
    hm = heatmap!(ax,qx,qy, Sq.(qs))
    Colorbar(fig[1,2],hm,label = L"\langle S^z(\mathbf{q})^*  S^z(\mathbf{q})\rangle")
    fig
end