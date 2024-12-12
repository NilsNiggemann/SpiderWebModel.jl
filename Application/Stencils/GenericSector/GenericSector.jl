cd(@__DIR__)
using Pkg
Pkg.activate("../..")
using CairoMakie, MakieHelpers,Statistics
import SpiderWebModel as SW

include("../plottingUtils.jl")
SW.initGurobi()
##

# SW.Random.seed!(1233)
SW.Random.seed!(100)
L = 20
@time sols = SW.constructGroundstatesSpin1(L, 50000, 0.35,boundaryCondition = :periodic,STotZero=true,progress = true,TimeLimit = 100)
##
confs = [SW.stencilConfig(float(x),1,boundaryCondition = :periodic) for x in sols.solutions]
##
sort!(confs,by = SW.NPlaquettes,rev=true)
S = confs[1]
SW.plotApplPlaquettes(S)
##
SW.Random.seed!(1234)
mu = 0.2

CTSR = SW.ContinuousTimeMethod(10,Hxx = SW.Hxx_RK(mu))


nThermal = 1000

ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([2,2],[2,-2]),S)
SW.rand!(ψGSymm,1e-4)
w_avg_estimate = 0.1*length(S)
##
SRSteps = 250
##
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,40 ,ψGSymm,SRSteps,8e-3,SW.IterativeSRSolver();Nwalkers = 3*20,reconfigure=false,rel_tolerance=0,equilibration_steps=nThermal,pre_equilibration_steps=40_000,
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
SqsGFMC = SW.getSqsGFMC(results,1:50:100)
##


with_theme(theme_SimpleTicks()) do 
    SqMat = dropmean(SqsGFMC,dims=4)[:,:,end,:]
    SqErr = dropstd(SqsGFMC,dims=4)[:,:,end,:]
    fittingCoefs = optimizeCoeffs(SqMat)
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

    # for (phi,color) in zip([0,pi/4],colors)
    #     qpoints_raw = q_path.(qr,phi)
    #     qpoints = sort!(unique!(roundToTrueMomenta.(qpoints_raw,size(SqMat,1)-1)), by = SW.norm)

    #     Sqcut = Sq.(qpoints)
    #     Sqerrcut = Sqerr.(qpoints)
        
    #     # SqFT = [SqFieldTheory(q,1,10) for q in qpoints]
    #     SqFT = [SqFieldTheory(q,fittingCoefs...) for q in qpoints]
    #     scatter!(ax,qpoints,marker = '×' ,color = color)
    #     scatterlines!(axFT,Point.(qpoints),color = color,linestyle = :dash,marker = '●',markersize = 4)
    #     qnorms_sq = SW.norm.(qpoints).^2
    #     scatter!(ax2,qnorms_sq,Sqcut,
    #     marker = '×',markersize = 15,color = color)
    #     errorbars!(ax2,qnorms_sq,Sqcut,Sqerrcut,color = color,whiskerwidth = 6,linewidth=0.5)
    #     scatterlines!(ax2,qnorms_sq,SqFT,color = color,linestyle = :dash,marker = '●',markersize = 4)
    # end
    rowsize!(fig.layout,1,Relative(0.5))

    # text!(axFT,Point(pi,1.4pi),text=L"r = %$(strd(fittingCoefs[2]))",color = :white,align = (:center,:center))
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-10,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-10,0))
    # Label(fig[2,1, TopLeft()],L"c)$$",padding = (-30,0,-10,0))
    # Label(fig[3,1, TopLeft()],L"d)$$",padding = (-30,0,-10,0))

    fig
end

##
using CairoMakie
thetas = pi.*rand(100000)
phis = 2pi.*rand(100000)

sphere(theta,phi) = SA[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]

points = sphere.(thetas,phis)

hist(getindex.(points,2),bins = 500)

##
function uniformOnSphere(rng = SW.Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng)
    z = 2.0 * rand(rng) - 1.0
    r = sqrt(1.0 - z * z)
    return SA[r * cos(phi), r * sin(phi), z]
end

##
function uniformOnSphere()
    phi = 2.0 * pi * rand()
    z = 2.0 * rand() - 1.0
    r = sqrt(1.0 - z * z)
    return SA[r * cos(phi), r * sin(phi), z]
end
##

@time pointsU = [uniformOnSphere() for _ in 1:200^2*2]
##
with_theme(theme_SimpleTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    hist!(ax,getindex.(pointsU,1),bins = 100,label = "x", color = (:red,0.4))
    hist!(ax,getindex.(pointsU,2),bins = 100,label = "y", color = (:green,0.4))
    hist!(ax,getindex.(pointsU,3),bins = 100,label = "z", color = (:blue,0.4))
    axislegend(ax)
    fig
end