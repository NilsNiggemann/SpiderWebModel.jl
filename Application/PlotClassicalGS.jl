import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
using HDF5
cd(@__DIR__)
##
files = readdir("../Application/ConfsRaw/GurobiConfs",join=true)

Configs = 
    let 
    Confs = Vector{Matrix{Bool}}()

    for file in files
        h5open(file,"r") do f
            for k in keys(f["Confs"])
                confs = read(f["Confs/$k"])
                for conf in eachslice(confs,dims=3)
                    push!(Confs, collect(conf))
                end
            end
        end
    end
    Spinconfs = SW.floatSpinConfig.(Confs,0.5)
    @assert all(SW.fulFillsConstraint,Spinconfs)
    Spinconfs
end
##

SW.plotFractons(Configs[10])
##
function SqLargeN(qx,qy)
    cx = cos(qx)
    cy = cos(qy)
    sx = sin(qx)
    sy = sin(qy)
    c2x = cos(2*qx)
    c2y = cos(2*qy)

    num = 2*(cx -cy +2sx*sy)^2
    denom = 4 - 4*cx*cy - c2x*(1-2c2y) - c2y
    return num/denom
end
##
Spin1Configs = let
    files = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/DataRK/",join=true)
    filter!(contains("L=80"),files)

    c1 = SW.SpinConfig.(eachslice(h5read(first(files),"confs") ./2,dims=3),1)
    for f in files[2:end]
        c2 = SW.SpinConfig.(eachslice(h5read(f,"confs") ./2,dims=3),1)
        append!(c1,c2)
    end
    c1
end
##
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
##
Sq = SW.getEqualWeightStructureFac(Configs)

resRK = SW.readResults(first(readdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataRK/L=100/2/",join=true)),1000);

SqSpin1 = SW.getSqGFMC(resRK[end],1) ./4

##
with_theme(theme_PiTicks()) do
    fig = Figure(size = 0.8 .*(750,300))
    ticks = PiTicks([0,pi])

    ax = Axis(fig[1, 1], aspect = 1,xticks = ticks,yticks = ticks,xminorticksvisible = true ,xlabel = L"q_x",ylabel = L"q_y",yminorticksvisible = true,title = L"Gaussian approximation$$")
    ax2 = Axis(fig[1, 2], aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,xticks = ticks,yticks = ticks,title = L"class. spin-$1/2$")
    ax3 = Axis(fig[1, 3], aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,xticks = ticks,yticks = ticks,title = L"class. spin-$1$")

    SqFunc = SW.getSqCont(real(Sq.Sq))

    Sq1Func = SW.getSqCont(real(SqSpin1) ./2) # normalize by S(S+1)
    
    kx = collect(trueMomenta(-0.5pi,1.5pi,size(Sq.Sq,1)-1))
    # kx = kx
    ky = kx
    kxSpin1 = collect(trueMomenta(-0.5pi,1.5pi,size(SqSpin1,1)-1))
    # kxSpin1 = 2pi .* (0.5:1)
    kySpin1 = kxSpin1

    SqLN = [SqLargeN(kx,ky)/4 for kx in kxSpin1 , ky in kySpin1]
    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]
    SqSpin1Plot = [Sq1Func(kx,ky) for kx in kxSpin1 , ky in kySpin1]
    @info "" sum(SqLN)/length(SqLN) sum(SqPlot)/length(SqPlot) sum(SqSpin1)/length(SqSpin1) length(kxSpin1)
    cRange = extrema(SqPlot)

    hm = heatmap!(ax,kx,ky,SqLN,colorrange = cRange)
    hm2 = heatmap!(ax2,kx,ky,SqPlot,colorrange = cRange)
    hm3 = heatmap!(ax3,kx,ky,SqSpin1Plot ,colorrange = cRange)
    # hm = halfhalfheatmap!(ax,kx,ky,SqFunc,SqLargeN,x->-x+3pi,normalize = true)
    Colorbar(fig[1, 4], hm2,ticks = [0,0.2,0.4,0.6],height = Relative(0.9),label = L"\mathcal{S}(q)/S(S+1)")
    save("../Application/figs/Sq_comparison.png", fig,px_per_unit=3)
    # Colorbar(fig[1, 2], hm)
    fig
    
end
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1)
    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqRK,1)-1)
    SqLN = [SqLargeN(kx,ky) for kx in kx , ky in ky]
    
    hm = heatmap!(kx,ky,SqLN)
    Colorbar(fig[1,2],hm)
    save("Application/TalkFigs/SqLN.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1)
    SqFunc = SW.getSqCont(real(Sq.Sq))

    kx =ky= collect(trueMomenta(-0.5pi,1.5pi,size(Sq.Sq,1)-1))

    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]

    # kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqRK,1)-1)
    hm = heatmap!(kx,ky,SqPlot)
    Colorbar(fig[1,2],hm)
    save("Application/TalkFigs/SqIsing.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1)
    Sq = SW.getSqCont(SqRK)
    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqRK,1)-1)
    hm = heatmap!(kx,ky,Sq.(Iterators.product(kx,ky)))
    Colorbar(fig[1,2],hm)
    save("Application/TalkFigs/SqRK.pdf",fig)
    fig
end

##
import SpiderWebModel as SW
using CairoMakie
using MakieHelpers
import SpiderWebModel: getStairCase, StaticArrays
##
function getPeriodic(parent)
    state = parent |> Array
    SW.SpinConfig(SW.PeriodicMatrix(state), parent.S)
end
##

function solveED(state, args...; kwargs...)
    @time HilbertSpace = SW.generateHilbertSpace(state)
    @time sol = SW.SolveHKrylov(HilbertSpace.H, args...; kwargs...)
    E0 = sol.values[1]
    ψ0 = sol.vectors[1]
    return (;E0,ψ0,HilbertSpace)
end

function getObservables(res, sol)
    # state = res.AllStates[1].parent

    AllStates = fetch.([Threads.@spawn SW.spinConfig(state) for state in res.AllStates])

    m = SW.getMagnetization(AllStates, sol)

    mavg = sum(abs, m) / length(m)
    mavgBulk = sum(abs, m[3:(end-2), 3:(end-2)]) / length(m[3:(end-2), 3:(end-2)])
    magnetization = (; m, mavg, mavgBulk)

    structureFac = SW.getStructureFac(AllStates, sol)

    return (; magnetization..., structureFac...)
end
##
function plotMagnetization!(ax, observables)
    heatmap!(ax, observables.m, colormap = :grays)
end

function plotStructureFac!(ax, observables)
    (; k, Sq_k) = observables
    heatmap!(ax, k, k, real(Sq_k))
end

function plotOverview(Observables; title = L"")
    with_theme(theme_SimpleTicks()) do
        fig = Figure(; size = 0.8 .* (450, 600))

        axmag = Axis(fig[1, 1]; title, xlabel = L"x", ylabel = L"y", aspect = 1)
        axSq = Axis(
            fig[2, 1],
            xlabel = L"q_x",
            ylabel = L"q_y",
            aspect = 1,
            xticks = PiTicks(0:(0.5pi):(2pi)),
            yticks = PiTicks(0:(0.5pi):(2pi)),
        )

        hm = plotMagnetization!(axmag, Observables)
        Colorbar(fig[1, 2], hm, label = L"\langle S^z \rangle")
        hm = plotStructureFac!(axSq, Observables)
        Colorbar(fig[2, 2], hm, label = L"\mathcal{S}^{zz}(\mathbf{q})")
        fig
    end
end

getRKWavefunction(ψ) = SW.normalize!(one.(ψ))
plotOverview(res, sol; kwargs...) = plotOverview(getObservables(res, sol); kwargs...)

##
InitialConf = getPeriodic(SW.getStairCase(14))
SolStair = solveED(InitialConf)
##
States = SW.spinConfig.(SolStair.HilbertSpace.AllStates,Ref(InitialConf),Ref(SolStair.HilbertSpace.plaqMapping))

SqED = SW.getStructureFac(States, SolStair.ψ0)
SqRK = SW.getStructureFac(States, getRKWavefunction(SolStair.ψ0))

##
with_theme(theme_PiTicks()) do
    fig = Figure(size = 0.7 .*(600,480))
    ticks = PiTicks([0,pi])

    subgl_top = GridLayout()
    subgl_bot = GridLayout()

    axkwargs = (;    xminorticksvisible=true,yminorticksvisible=true,xtickwidth = 1.4,ytickwidth = 1.4,xminortickwidth = 1.2,yminortickwidth = 1.2,xticksize=5,yticksize=5,xminorticksize=3,yminorticksize=3,xminortickalign =1,yminortickalign =1,xtickalign =1,xticksmirrored =true,yticksmirrored=true,ytickalign=1,
    xlabelpadding = -5
    )

    subgl_top[1, 1] = ax = Axis(fig, aspect = 1,xticks = ticks,yticks = ticks,xminorticksvisible = true ,xlabel = L"q_x",ylabel = L"q_y",yminorticksvisible = true,
    xlabelvisible=true;
    axkwargs...
    )
    subgl_top[1, 2] = ax2 = Axis(fig, aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,
    xlabelvisible=true,
    xticklabelsvisible=true,xticks = ticks,yticks = ticks,
    # title = L"class. spin-$1/2$"
    ;
    axkwargs...
    )
    subgl_top[1, 3] = ax3 = Axis(fig, aspect = 1,yticklabelsvisible=false,xminorticksvisible = true ,xlabel = L"q_x",yminorticksvisible = true,
    xlabelvisible=true,
    xticklabelsvisible=true,xticks = ticks,yticks = ticks,
    # title = L"class. spin-$1$"
    ;
    axkwargs...
    )

    SqFunc = SW.getSqCont(real(Sq.Sq) ./ (1/2*(1/2+1)))

    Sq1Func = SW.getSqCont(real(SqSpin1) ./2) # normalize by S(S+1)
    
    kx = collect(trueMomenta(-0.5pi,1.5pi,size(Sq.Sq,1)-1))
    # kx = kx
    ky = kx
    kxSpin1 = collect(trueMomenta(-0.5pi,1.5pi,size(SqSpin1,1)-1))
    # kxSpin1 = 2pi .* (0.5:1)
    kySpin1 = kxSpin1

    SqLN = [SqLargeN(kx,ky)/4 for kx in kxSpin1 , ky in kySpin1]
    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]
    SqSpin1Plot = [Sq1Func(kx,ky) for kx in kxSpin1 , ky in kySpin1]

    @info "" sum(SqLN)/length(SqLN) sum(SqPlot)/length(SqPlot) sum(SqSpin1)/length(SqSpin1) length(kxSpin1)
    # cRange = extrema(SqPlot)
    cRange = (0,maximum(SqPlot))
    hm = heatmap!(ax,kx,ky,SqLN,colorrange = cRange)
    hm2 = heatmap!(ax2,kx,ky,SqPlot,colorrange = cRange)
    hm3 = heatmap!(ax3,kx,ky,SqSpin1Plot,colorrange = cRange)
    # hm = halfhalfheatmap!(ax,kx,ky,SqFunc,SqLargeN,x->-x+3pi,normalize = true)

    subgl_bot[1, 1] = axED1 = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
    # axkwargs...,
    # xtickcolor = :white,ytickcolor = :white,xminortickcolor = :white,yminortickcolor = :white
    )
    subgl_bot[1, 2] = axED2 = Axis(fig, xlabel = L"q_x", ylabel = L"q_y", aspect = 1,ylabelvisible = false,yticklabelsvisible = false,xminorticksvisible = true, yminorticksvisible = true,xticks = ticks, yticks = ticks;
    # axkwargs...,
    # xtickcolor = :white,ytickcolor = :white,xminortickcolor = :white,yminortickcolor = :white
    )
    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqED.Sq,1)-1)
    SqFunc = SW.getSqCont(SqED.Sq  ./ (1/2*(1/2+1)))
    SqMat = [real(SqFunc(x,y)) for x in kx, y in ky] 

    colorrange = extrema(SqMat)
    SqRKFunc = SW.getSqCont(SqRK.Sq  ./ (1/2*(1/2+1)))
    SqRKMat = [real(SqRKFunc(x,y)) for x in kx, y in ky]
    hmED1 = heatmap!(axED1, kx, ky, SqMat;colorrange )
    colorrange = extrema(SqRKMat)
    hmED2 = heatmap!(axED2, kx, ky, SqRKMat;colorrange)

    fig.layout[1, 1] = subgl_top
    fig.layout[2, 1] = subgl_bot
    # Label(fig[1,1, TopLeft()],L"a)$$",padding = (-30,0,-20,0))
    # Label(fig[1,2, TopLeft()],L"b)$$",padding = (-30,0,-20,0))
    # Label(fig[1,3, TopLeft()],L"c)$$",padding = (-30,0,-20,0))
    # Label(fig[2,1, TopLeft()],L"d)$$",padding = (-30,0,-20,0))
    # Label(fig[2,2, TopLeft()],L"e)$$",padding = (-30,0,-20,0))

    textpos = Point(-pi/2,3pi/2)


    Colorbar(subgl_top[1,4], hm2,height = Relative(1),width = Relative(0.8),label = L"\mathcal{S}(q)/S(S+1)",ticks = SimpleTicks([0,0.2,0.4,0.6]))
    Colorbar(subgl_bot[1,3], hmED2,height = Relative(1),width = Relative(0.8),label = L"\mathcal{S}(q)/S(S+1)",ticks = SimpleTicks())


    text!(ax, textpos ,text = L"a)",color = :black,align = (:left,:top),fontsize = 18)
    text!(ax2, textpos ,text = L"b)",color = :black,align = (:left,:top),fontsize = 18)
    text!(ax3, textpos ,text = L"c)",color = :black,align = (:left,:top),fontsize = 18)

    text!(axED1, textpos ,text = L"d)",color = :white,align = (:left,:top),fontsize = 18)
    text!(axED2, textpos ,text = L"e)",color = :white,align = (:left,:top),fontsize = 18)

    rowsize!(fig.layout,2,Relative(0.6))
    rowgap!(fig.layout,1,-5)

    colsize!(subgl_top,4,Relative(0.02))
    colsize!(subgl_bot,3,Relative(0.05))

    colgap!(subgl_top,1,0)
    colgap!(subgl_top,2,0)
    colgap!(subgl_top,3,5)
    save("../Application/figs/Sq_comparison_2.png", fig,px_per_unit=4)
    # Colorbar(fig[1, 2], hm)
    fig
    
end