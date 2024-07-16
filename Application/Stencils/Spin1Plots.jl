using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
using Optim
cd(@__DIR__)
function trueMomenta(kmin,kmax,L)
    nmin = floor(Int,L*kmin/(2pi))
    nmax = ceil(Int,L*kmax/(2pi))
    # return 1/(2pi*L*100) .* nmin:nmax
    return (nmin : nmax) .* 2pi/L
end
##
# filesRK = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/DataRK/",join=true)
# resRK = h5open(filesRK[end]) do f
#     read(f)
# end
##
resRK = SW.readResults(first(readdir("/scratch/hpc-prf-pm2frg/niggeni/Spiderweb/DataRK/L=100/2/",join=true)),1000);
##
SqRK = let  
    # RKConfs = eachslice(resRK[1].SaveConfigs,dims=(3,4))
    SqRK = SW.getSqGFMC(resRK[end],1)
    
end

##
res = h5open("../Data/Spin1GFMC_Eval_periodic.h5") do f
    read(f)
end

resOld = h5open("../Data/Spin1GFMC_Eval_periodic_L40.h5") do f
    read(f)
end

projectionSteps = 500
SqsGFMC = eachslice(resOld["SqsGFMC"][string(projectionSteps)],dims=3)

SqMat = mean(real(SqsGFMC)) ./4
errMat = std(real(SqsGFMC)) ./4
SqFunc = SW.getSqCont(SqMat)

errFunc = SW.getSqCont(errMat)

##
with_theme(theme_SimpleTicks()) do
    L = 20
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0/L^2",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    
    for (L,resL) in res
        L = parse(Int,L)
        en = resL["energies"]
        proj = resL["nBra"] .*axes(en,1)
        en = en[1:2:24,:]
        proj = proj[1:2:24]
        e0s,e0s_err = mean(en,dims=2)[:],sqrt.(var(en,dims=2)[:])
        e0s ./= L^2
        e0s_err ./= L^2
        # e0s = (e0s .- minimum(e0s) .+1e-1)
        
        ll = scatterlines!(ax,proj,e0s,label = L"GFMC $L = %$L$", marker = '●',markersize = 5)
        # errorbars!(ax,proj,e0s,e0s_err,whiskerwidth = 3.5)
        col = ll.color[]
        band!(ax,proj,e0s .- e0s_err,e0s .+ e0s_err,color = (col,0.2))
    end
    axislegend(ax,merge=true)
    ylims!(ax,-0.266,-0.265)
    # hlines!(ax,[-0.26582],color=:black,linestyle = :dash)
    # save("../figs/PaperFigs/GFMCEnergy.png",fig)
    fig
end
##
function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end
SqFieldTheory(k) = SqFieldTheory(k[1],k[2])
##
function SqFieldTheory2(kx,ky,k,b2,b3,A=1)
    A*(sqrt(2)*sqrt(-((-4 + 4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky))*(40*(b2 + b3) + k + 8*b2*(-8*cos(kx)*cos(ky) + cos(2*ky)) + 8*cos(2*kx)*(b2 - 4*b3 + (b2 + 2*b3)*cos(2*ky)) + 4*b3*(cos(4*kx) - 8*cos(2*ky) + cos(4*ky)))))*(cos(kx) - cos(ky) + 2*sin(kx)*sin(ky))^2)/((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2+1e-30)
end

SqFieldTheory2(k::AbstractVector,a,b,c,A=1) = SqFieldTheory2(k[1],k[2],a,b,c,A)

function SqFieldTheory3(kx,ky,b2,b3)
    SqFieldTheory2(kx,ky,0,b2,b3,1)
end

SqFieldTheory3(k::AbstractVector,b,c) = SqFieldTheory3(k[1],k[2],b,c)
##
function optimizeCoeffs(SqMat)
    k = 2pi .* LinRange(0,1,size(SqMat,1))

    # SqFT = [SqFieldTheory2(x,y,0,0) for x in kx, y in ky]
    function loss(b2,b3)
        l = 0.
        b2 = abs(b2)
        b3 = abs(b3)
        for (i,kx) in enumerate(k), (j,ky) in enumerate(k)
        # for (i,kx) in enumerate(k)
            # ky = 0
            # j = 1
            # l += abs2(SqMat[i,j] - SqFieldTheory2(kx,ky,v,w))#/(SqMat[i,j]+1e-5)
            # if (kx^2 +ky^2) < (0.8pi)^2
            # l += abs2(SqMat[i,j] - SqFieldTheory2(kx,ky,b1,b2,b3))*(1/(SqMat[i,j]+1e-10) + 1)
            l += abs2(SqMat[i,j] - SqFieldTheory3(kx,ky,b2,b3))#*(0.7/(kx^2 +ky^2 +1e-2) + 1)
            # end
        end
        # for (i,kx) in enumerate(kx)
        #     l += abs2(SqMat[i,1] - SqFieldTheory2(kx,0,v,w))/(SqMat[i,1]+1e-5)
        # end
        return l
    end
    loss(v) = loss(v[1],v[2])

    x0 = [1, 1.]

    res = optimize(loss, x0)
    @info res
    coefs = abs.(Optim.minimizer(res))
    return coefs
end

fittingCoefs = optimizeCoeffs(SqMat)
##
function getFittedSq(coefs)
    SqFit(kx,ky) = SqFieldTheory3(kx,ky,coefs...)
    SqFit(k::Union{AbstractVector,Tuple}) = SqFieldTheory3(k[1],k[2],coefs...)
    SqFit
end
SqFieldTheoryFit = getFittedSq(fittingCoefs)
# fittingCoeffs::NamedTuple{(:v, :w), Tuple{Float64, Float64}} = optimizeCoeffs(SqMat)
##

# SqFieldTheory(k) = SqFieldTheory2(k,1/850,0.2/850)
# SqFieldTheory(kx,ky) = SqFieldTheory2(kx,ky,1/850,0.2/850)


with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx =ky= collect(trueMomenta(-0.5pi,1.5pi,40))
    
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1)

    SqFTfunc = getFittedSq(fittingCoefs)
    # SqFTfunc = getFittedSq([1,0.1,0.1])
    SqFT = [SqFTfunc(x,y) for x in kx, y in ky]
    hmFT = heatmap!(ax,kx,ky,SqFT,colormap = :viridis)
    
    Colorbar(fig[1,2],hmFT,
    # label = L"\mathcal{S}(\textbf{q})",
    height = Relative(0.95),vertical=true,ticks = SimpleTicks(),flipaxis=true)
    # save("../figs/TalkFigs/SqFT_b2_b3_fit.pdf",fig)
    # save("../figs/TalkFigs/SqFT_B_0_0.pdf",fig)
    fig
end
##

with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 20,size = 2 .*(450,155))
    # fig = Figure()
    
    ticks = PiTicks([0,pi,])
    Ls = sort(collect(keys(res)))

    minval = Inf
    maxval = -Inf
    projectionSteps = 500
    SqFTplot = getFittedSq(fittingCoefs)
    for (L,resL) in res
        SqsGFMC_L = eachslice(resL["SqsGFMC"][string(projectionSteps)],dims=3)

        SqMat = mean(real(SqsGFMC_L)) ./4

        minval = min(minval,minimum(SqMat))
        maxval = max(maxval,maximum(SqMat))
    end

    for (i,Lstr) in enumerate(Ls)
        resL = res[Lstr]
        L = parse(Int,Lstr)
        SqsGFMC_L = eachslice(resL["SqsGFMC"][string(projectionSteps)],dims=3)

        SqMat = mean(real(SqsGFMC_L)) ./4
        SqMat ./= maxval

        # errMat = sqrt.(var(real(SqsGFMC))) ./4
        SqFunc = SW.getSqCont(SqMat)

        kx = ky = trueMomenta(-0.5pi,1.5pi,L)
        Sq = [SqFunc(x,y) for x in kx, y in ky]

        axFT = Axis(fig[1,i],xlabel = L"q_x",ylabel = L"q_y",title = L"$L=%$L$",aspect = 1,xticks = ticks,yticks = ticks,xlabelvisible=true,yticklabelsvisible=i==1,ylabelvisible = i==1,xticklabelsvisible=true)

        hmMC = halfhalfheatmap!(axFT,kx,ky,SqFTplot,SqFunc,x->-x+pi,normalize = true)

    end
    Colorbar(fig[:,end+1],colorrange =(minval,maxval),)
    # rowgap!(fig.layout,1,-6)    
    colgap!(fig.layout,1,-4)
    colgap!(fig.layout,2,-4)
    colgap!(fig.layout,3,-4)
    text!(fig[1,1],(-0.4pi,-0.4pi),text = L"U(1)$$",color = :white,align = (:left,:bottom))
    text!(fig[1,1],(1.4pi,1.4pi),text = L"$S=1$",color = :white,align = (:right,:top))
    # colgap!(fig.layout,4,-4)
    # rowsize!(fig.layout,2,Relative(0.05))
    save("../figs/PaperFigs/Sq_size_scaling.pdf",fig)
    fig
end
##
with_theme(theme_SimpleTicks()) do 
    kx = trueMomenta(-0.5pi,1.5pi,30)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"\mathcal{S}(\mathbf{q})",xticks = PiTicks())
    
    colors = (:black,:lightblue,:darkblue,)
    for (color,p) in zip(colors,(100,200,500))
        SqsGFMC_p = eachslice(res["30"]["SqsGFMC"][string(p)],dims=3)

        SqMat = mean(real(SqsGFMC_p)) ./4
        errMat = std(real(SqsGFMC_p)) ./4
        SqFunc = SW.getSqCont(SqMat)

        errFunc = SW.getSqCont(errMat)

        Sqerr = [errFunc(x,0) for x in kx]
        SqVec = [SqFunc(x,0) for x in kx]
        # SqVec = SqMat[1,:]
        # Sqerr = errMat[1,:]
        # SqVec ./= maximum(SqVec)
        # Sqerr ./= maximum(SqVec)

        scatterlines!(kx,SqVec;color)
        errorbars!(kx,SqVec,Sqerr;color,whiskerwidth = 8)
    end
    fittingCoefs2 = (fittingCoefs[1:end-1]./350000...,1) #./ last(fittingCoefs) #./630 
    SqFTfunc = getFittedSq(fittingCoefs)
    SqFT = [SqFTfunc(x,0) for x in kx]
    # SqFT ./= maximum(SqFT)

    lines!(kx,SqFT,color = :red)
    fig
end
##

with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    # kxSmall = kySmall = filter(x->abs(x) < pi/2,kx)
    kxSmall = kySmall = kx

    err = [errFunc(x,y) for x in kxSmall, y in kySmall]
    fig = Figure(fontsize = 20,size = (600,400))
    
    ticks = PiTicks([0,pi])
    # ticks = ([-3pi/2,-pi/2,pi/2],[L"-3π/2",L"-π/2",L"π/2"])
    # smallticks = ([-pi/3,0,pi/3],[L"-π/3",L"0",L"π/3"])
    smallticks = ticks

    axMC = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,xticks = ticks,yticks = ticks)
    
    
    axerr = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=true,xticks = smallticks,yticks = smallticks)

    hmMC = halfhalfheatmap!(axMC,kx,ky,SqFieldTheory,SqFunc,x->-x-pi,normalize = true)
    # SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    # heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(!isnan,Sq))
    hmerr = heatmap!(axerr,kxSmall,kySmall,err,colormap = :viridis)
    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)

    Colorbar(fig[2,1],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,2],hmerr,label = L"\sigma\mathcal{S}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    
    # colsize!(fig.layout,1,Relative(0.6))
    rowsize!(fig.layout,2,Relative(0.1))
    fig
end
##
function rasterCurve(curvePoints,grid,t)

    getPos(point) = findmin(x->SW.norm(SW.SVector(x.-point)),grid)[2]
    positions = getPos.(curvePoints)
    tnew = empty(t)
    posnew = empty(positions)
    for i in eachindex(t)
        p = positions[i]
        if p ∉ posnew
            push!(tnew, t[i])
            push!(posnew,p)
        end 
    end
    return tnew,posnew
end

using StaticArrays
function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,SW.norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end
##
with_theme(theme_PiTicks()) do 

    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = (650,450))
    
    ticks = PiTicks([0,pi,])

    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    axMC = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    # axerr = Axis(fig[1,3],xlabel = L"q_x",ylabel = L"q_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    
    axpath = Axis(fig[2,1:3],xlabel = L"t",ylabel = L"\tilde{\mathcal{S}}(\mathbf{q})",yticks = SimpleTicks(),yminorticksvisible = true,xminorticksvisible = true)

    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheoryFit(x,y) for x in kx, y in ky]

    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))

    tRange = LinRange(0,2pi,100)
    path(r,t) = (r*cos(t),r*sin(t))

    p1 = [path(0.5,t) for t in tRange]
    xygrid = [(x,y) for x in kx, y in ky]

    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    # filter!(x -> x[1] in kx && x[2] in ky,p1)
    lines!(axFT,Point.(xygrid[p1_discrete]),color = :red,linewidth = 2)
    lines!(axMC,Point.(xygrid[p1_discrete]),color = :black,linewidth = 2)
    Sqp = Sq[p1_discrete]
    
    Sqp ./= maximum(Sqp)
    # SqFT_coarse = SW.getSqCont(SqFT)
    Sqp_all = [S[p1_discrete]./4 for S in real(SqsGFMC)]
    Sqerr = sqrt.(var(Sqp_all))
    # SqFTp = SqFieldTheory.(p1)
    SqFTp = SqFT[p1_discrete]
    SqFTp ./= maximum(SqFTp)
    scatterlines!(axpath,tRange,SqFTp,color = :red,linewidth = 1,linestyle = :solid)
    scatter!(axpath,tRange,Sqp,color = :black,markersize = 8)
    errorbars!(axpath,tRange,Sqp,Sqerr,color = :black,whiskerwidth = 8)
    Colorbar(fig[1,3],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(1),ticks = SimpleTicks())
    rowgap!(fig.layout,1,Relative(-0.))
    rowsize!(fig.layout,2,Relative(0.5))
    xlims!(axpath,0,2pi)
    # save("../figs/PaperFigs/Spin1/Sq_L40.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do 

    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    # Sq ./= maximum(Sq)
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = 1.1 .*(480,450))
    
    ticks = PiTicks([0,pi,])

    axMC = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",title = L"spin model$$",aspect = 1,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -10,
    ylabelvisible = false,yticklabelsvisible=false,
    )
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"U(1) theory$$",aspect = 1,
    xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -10,
    )
    
    axpath1 = Axis(fig[2,1:3],xlabel = L"k^4",ylabel = L"\mathcal{S}(\mathbf{q})",yticks = SimpleTicks(),xticks = PiTicks(),yminorticksvisible = true,xminorticksvisible = true,xminorgridvisible =true,yminorgridvisible=true,yaxisposition=:left,

    )
    axpath2 = Axis(fig[2,1:3],
    yaxisposition=:right,yticklabelcolor=:red,
    yticks = SimpleTicks(),
    xlabelvisible = false,
    ygridvisible=false,
    xgridvisible=false,
    xticklabelsvisible= false,
    xticksvisible= false,
    # xminorticksvisible=false
    ylabelvisible=false,
    )

    # SqFTPl(x,y) = SqFieldTheoryFit(x,y)
    SqFTPl = getFittedSq(fittingCoefs)

    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFTPl(x,y) for x in kx, y in ky]
    # SqFT ./= maximum(SqFT)

    hmFT = heatmap!(axFT,kx,ky,SqFT ,colorrange = extrema(Sq),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    path(t,angle) = t .*(cos(angle),sin(angle)) #.+ (pi,pi)
    
    
    tRange = LinRange(0,.3pi,500)
    # angle = deg2rad(45)
    angle = pi/3
    p1 = [path(t,angle) for t in tRange]
    xygrid = [(x+0.5pi,y+0.5pi) for x in kx, y in ky]
    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    tRange² = (tRange) .^4
    # filter!(x -> x[1] in kx && x[2] in ky,p1)
    lines!(axFT,Point.(xygrid[p1_discrete]);color=:red,linewidth = 2,linestyle = :dash)
    lines!(axMC,Point.(xygrid[p1_discrete]);color=:black,linewidth = 2)
    colors = (:black,:lightblue,:darkblue)
    for (color,p) in zip(colors,(500,200))
        SqsGFMC = eachslice(res["30"]["SqsGFMC"][string(p)],dims=3)./4

        Sqp_all = [S[p1_discrete] for S in real(SqsGFMC)]
    
        Sqp = mean(Sqp_all)
        Sqerr = sqrt.(var(Sqp_all))

        scatter!(axpath1,tRange²,Sqp;color,markersize = 8,label = "p = $p")
        errorbars!(axpath1,tRange²,Sqp,Sqerr;color,whiskerwidth = 8)
        # xlims!(axpath1,-0.1,0.85pi)
    end
    axislegend(axpath1,merge=true,position=:lt)
    SqFTp = [SqFTPl(x,y) for (x,y) in (xygrid[p1_discrete])]

    lines!(axpath2,tRange²,SqFTp;color=:red,linewidth = 2,linestyle = :dash)
    # lines!(axpath1,tRange²,tRange.^4*0.04/(0.03pi);color=:blue,linewidth = 2,linestyle = :solid)
    
    Colorbar(fig[1,3],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(1),ticks = SimpleTicks())
    Label(fig[1,1, TopLeft()],L"a)$$",padding =(-45,0,0,0))
    Label(fig[1,2, TopLeft()],L"b)$$",padding =(-20,0,0,0))
    Label(fig[2,1, TopLeft()],L"c)$$",padding =(-45,10,10,0))
    # rowgap!(fig.layout,1,Relative(-0.))
    rowgap!(fig.layout,1,Relative(-0.1))
    rowsize!(fig.layout,2,Relative(0.5))

    save("../figs/PaperFigs/Sq_L40_comp.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1)
    res40 = h5open("../Data/Spin1GFMC_Eval_open.h5") do f
        read(f)
    end



    # SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]

    SqsGFMC = mean(res40["SqsGFMC"]["750"],dims=3)./4

    kx =ky= collect(trueMomenta(-0.5pi,1.5pi,size(SqsGFMC,1)-1))

    SqFunc = SW.getSqCont(real(SqsGFMC))

    SqPlot = [SqFunc(kx,ky) for kx in kx , ky in kx]

    # kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqRK,1)-1)
    hm = heatmap!(kx,ky,SqPlot)
    Colorbar(fig[1,2],hm,ticks = SimpleTicks())
    save("../figs/TalkFigs/SqL40.pdf",fig)
    fig
end


##
function omegaFT(x,y)
    a = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return sqrt(a)
end
omegaFT(k) = omegaFT(k[1],k[2])


with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = 2pi .* LinRange(0,1,500) .- 0.5pi
    
    fig = Figure(fontsize = 20,size = (450,400))
    
    ticks = PiTicks([0,pi,])
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks)

    SqFT = [omegaFT(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    
    Colorbar(fig[1,2],hmFT,label = L"\omega(\textbf{q})",height = Relative(0.95),vertical=true,ticks = SimpleTicks(),flipaxis=true)
    save("../figs/PaperFigs/omegaq_fieldTheory.png",fig)
    fig
end
##
BBCorrsFile = "../Data/BBCorr2/BBCorr_L=20.h5"
BBCorrs = h5read(BBCorrsFile,"BBCorrelator")
L = h5read(BBCorrsFile,"L")
S = SW.stencilConfig(zeros(L,L),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())

nBra = only(unique!(h5read(BBCorrsFile,"nBra")))

BBCorrelator_p = mean(BBCorrs,dims=1)[1,:,:]
BBCorrelator_err_p = sqrt.(var(BBCorrs,dims=1))[1,:,:]
refPlaq = SW.getCentralPlaquette(S)
##
AllPlaqs = collect(SW.plaquetteIterator(S))
AllPlaqsDict = Dict(AllPlaqs[i] => i for i in eachindex(AllPlaqs))

with_theme(theme_SimpleTicks()) do 
    fig = Figure(fontsize = 20)
    ax = Axis(fig[1,1],xlabel = L"\mathcal{P}",ylabel = L"\mathcal{B^2}_{i,j}",yticks = SimpleTicks(),yminorticksvisible = true,xminorticksvisible = true)
    
    projaxis = nBra .* axes(BBCorrelator_p,2)
    for P in ((0,0),(1,1),(2,0),(2,2),(3,1),(3,3))
    # for P in ((1,1),(2,0),(2,2),(3,1),(3,3),(4,0),(4,2),(4,4))

        i_plaq = AllPlaqsDict[refPlaq .+ P]

        lines!(projaxis,BBCorrelator_p[i_plaq,:])
        errorbars!(projaxis,BBCorrelator_p[i_plaq,:],BBCorrelator_err_p[i_plaq,:],whiskerwidth = 8)
    end
    fig    
end
##
BBCorrelator = BBCorrelator_p[:,5]
BBCorrelator_err = BBCorrelator_err_p[:,5]
# BBCorrelator = BBCorrs[10,1,:]
S = SW.stencilConfig(zeros(L,L),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())

##
# SW.plotApplPlaquettes(parentState); scatter!(Point.(GFMCPlaqs),markersize = 100 .* BBCorrelator); current_figure()
function plotRealSpaceCorr(S,plaqs,BBCorr,BBCorr_err)
    fig = SW.plotApplPlaquettes(S,color = (:red,0.0),markersize = 1)
    maxVal = maximum(abs,BBCorr)
    sizefunc(x) = min(2 * abs(x)/maxVal,50)

    function plotcircle!(r0,radius;kwargs...)
        r(x) = r0 .+ radius .* [cos(x),sin(x)]
        lines!(r.(LinRange(0,2pi,100));kwargs...)
    end
    scatter!(Point(refPlaq),marker = '×',color = :red,markersize = 25)
    # scatter!(Point.(plaqs),markersize = sizefunc.(BBCorr .- BBCorr_err),marker = '○')
    for (i,p) in enumerate(plaqs)
        l = plotcircle!(SVector(p),sizefunc(BBCorr[i] - BBCorr_err[i]);linewidth = 0.4)
        plotcircle!(SVector(p),sizefunc(BBCorr[i] + BBCorr_err[i]);linewidth = 0.4,color = l.color[])
    end
    # scatter!(Point.(plaqs),markersize = sizefunc.(BBCorr .- BBCorr_err),marker = '○')
    # scatter!(Point.(plaqs),markersize = sizefunc.(BBCorr .+ BBCorr_err),marker = '○')
    save("../figs/PaperFigs/Spin1/RealSpaceCorr_L20.pdf",fig)
    fig
end
plotRealSpaceCorr(S,AllPlaqs,BBCorrelator,BBCorrelator_err)
##
function FTPlaq(rPlaq,Vals,k)
    res = 0.
    for (r,Val) in zip(rPlaq,Vals)
        res += Val * cos(k'*r)
    end
    res
end

function ω_photon(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
ω_photon((kx,ky)) = ω_photon(kx,ky)

with_theme(theme_PiTicks()) do
    T = SW.SA[
        1 1;
        -1 1
    ]/2
    # ri = [T * SW.SVector(r .- refPlaq) for r in AllPlaqs]
    ri = [SW.SVector(r .- refPlaq) for r in AllPlaqs]
    # rPlaq = [mapToPlaquetteBasis(r) for r in ri]
    capfilter(x) = min(abs(x),30)


    # return scatter(Point.(ri),color = capfilter.(200 .*BBCorrelator./maximum(BBCorrelator)),markersize = capfilter.(400 * BBCorrelator./maximum(BBCorrelator)),axis =(;aspect=1,xticks=SimpleTicks(-L:L), yticks=SimpleTicks(-L:L)))

    # return scatter(Point.(ri),color = BBCorrelator,markersize = 150 * abs.(BBCorrelator))
    k = LinRange(-pi,pi,200)

    # FT = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    FTgen(BBCorrelator) = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    hm = nothing
    fig = Figure(size = 1 .*(500,600))

    for (i,p) in enumerate((1,5,8))
        FTs = FTgen.(eachslice(BBCorrs[:,:,p],dims=1))
        FT = mean(FTs)
        FTerr = sqrt.(var(FTs))#

        title1 = i == 1 ? L"err$$" : ""
        ax = Axis(fig[i,1];aspect = 1,ylabel = 
        
        L"\mathcal{P} = %$(nBra*p)")
        ax2 = Axis(fig[i,2];aspect = 1,title = title1)
        # ax3 = Axis(fig[1,3];aspect = 1)
        
        hm = heatmap!(ax,k,k,FT)
        hm2 = heatmap!(ax2,k,k,FTerr,colorrange = extrema(FT))
    end
    # photw = ω_photon.(Iterators.product(k,k))
    # photw .*= maximum(FT)/maximum(photw)

    # heatmap!(ax3,k,k,photw,colormap = :viridis,colorrange = extrema(FT))
    Colorbar(fig[:,3],hm,label = L"\mathcal{B}^2(\mathbf{q})",ticks = SimpleTicks())
    save("../figs/PaperFigs/Spin1/BCorrFT_L20_comp.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do 
    fig = Figure(size = 1 .*(500,600))
    k = LinRange(-pi,pi,200)
    photw = ω_photon.(Iterators.product(k,k))
    # photw .*= maximum(FT)/maximum(photw)
    ax = Axis(fig[1,1];aspect = 1)

    hm = heatmap!(ax,k,k,photw,colormap = :viridis)
    fig
end

##
BBMat,BBerr = let 
    T = SW.SA[
        1 1;
        -1 1
    ]/2
    ri = [ SW.SVector(r .- refPlaq) for r in AllPlaqs]
    # ri = [T * SW.SVector(r .- refPlaq) for r in AllPlaqs]
    # rPlaq = [mapToPlaquetteBasis(r) for r in ri]


    # k = LinRange(-pi,pi,200)
    k = trueMomenta(-pi,pi,20)

    FTgen(BBCorrelator) = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    hm = nothing
    p = 5
    FTs = FTgen.(eachslice(BBCorrs[:,:,p],dims=1))
    FT = mean(FTs)
    FTerr = std(FTs)
    FT,FTerr
end
# BBMat,BBerr = let 
#     ri = [ SW.SVector(r) for r in AllPlaqs]

#     FTMat = zeros(maximum(first,ri),maximum(last,ri))


#     function FTgen(BBCorrelator)
#         for (r,B) in zip(ri,BBCorrelator)
#             I = CartesianIndex(r[1],r[2])
#             FTMat[I] = B
#         end
#         return real(SW.FFTW.fft(FTMat))
#     end
#     p = 5
#     FTs = FTgen.(eachslice(BBCorrs[:,:,p],dims=1))
#     FT = mean(FTs)
#     FTerr = std(FTs)
#     FT,FTerr
# end

BBfunc = SW.getSqCont(BBMat)
BBerrfunc = SW.getSqCont(BBerr)
##
with_theme(theme_PiTicks()) do
    
    capfilter(x) = min(abs(x),30)


    k = LinRange(-pi,pi,200)

    hm = nothing
    fig = Figure(size = 1. .*(400,320))
    p = 5

    ax = Axis(fig[1,1];aspect = 1,ylabel = L"q_y",xlabel = L"q_x")
    # ax3 = Axis(fig[1,3];aspect = 1)
    
    hm = heatmap!(ax,k,k,BBMat)

    Colorbar(fig[1,2],hm,label = L"\mathcal{B}^2(\mathbf{q})",ticks = SimpleTicks())
    # save("../figs/PaperFigs/Spin1/BCorrFT_L20.pdf",fig)
    fig
end
##

with_theme(theme_PiTicks()) do 

    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    # Sq ./= maximum(Sq)
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = 1.1 .*(580,450))
    
    fig_top = fig[1,1] = GridLayout()
    fig_bottom = fig[2,1] = GridLayout()
    fig_bottom_left = fig_bottom[1,1:4] = GridLayout()
    fig_bottom_right = fig_bottom[1,5:6] = GridLayout()

    ticks = PiTicks([0,pi,])

    # kpointlabels = [L"Γ",L"X",L"X'",L"Γ"]
    KPoints = Dict([
        "Γ" => SVector(0,0),
        "X" => SVector(pi,0),
        "M" => SVector(pi,pi),
        "X'" => SVector(0,pi)
    ])
    kpath = ["Γ","X","X'","Γ"]
    pointlabels,p1 = fetchKPath([KPoints[k] for k in kpath],500)
    kpointlabels = Makie.latexstring.(kpath)

    # pointlabels,p1 = fetchKPath([SA[0,0],SA[pi,0],SA[pi,pi],SA[0.0,0.0]],500)
    # pointlabels,p1 = fetchKPath([SA[0,0.],SA[1,1.]],500)
    tRange = eachindex(p1)

    axFT = Axis(fig_top[1,1],xlabel = L"q_x",ylabel = L"q_y",
    # title = L"U(1) theory$$",
    aspect = 1,
    xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -10,
    )    

    axFTB = Axis(fig_top[1,2],xlabel = L"q_x",ylabel = L"q_y",
    # title = L"U(1) theory$$",
    aspect = 1,
    xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -10,
    ylabelvisible = false,yticklabelsvisible=false,
    )   
    
    axMC = Axis(fig_top[1,3],xlabel = L"q_x",ylabel = L"q_y",
    # title = L"spin model$$",
    aspect = 1,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -10,
    ylabelvisible = false,yticklabelsvisible=false,
    )
    axBB = Axis(fig_bottom_right[1,1],xlabel = L"q_x",ylabel = L"q_y",aspect = 1,
    xticks = PiTicks([-pi,0,pi]),yticks = PiTicks([-pi,0,pi]),yminorticksvisible = true,xminorticksvisible = true,
    xlabelpadding = -5,
    ylabelpadding = -10,
    title = L"\mathcal{B}^2(\mathbf{q})",
    # title = L"\langle \mathcal{F}(\mathbf{q}) \mathcal{F}^\dagger(-\mathbf{q})\rangle",
    )
    
    axpath1 = Axis(fig_bottom_left[:,1:3],ylabel = L"\mathcal{S}(\mathbf{q})",yticks = SimpleTicks(),
    # xticks = (tRange[pointlabels],[L"Γ",L"X",L"M",L"Γ"],)
    xticks = (tRange[pointlabels],kpointlabels,),
    # xticks = (tRange[pointlabels],[L"Γ",L"X"],),
    yminorticksvisible = true,xminorticksvisible = true,xminorgridvisible =true,yminorgridvisible=true,yaxisposition=:left,

    )
    # axpath2 = Axis(fig[2,2:4],
    # yaxisposition=:right,yticklabelcolor=:red,
    # yticks = SimpleTicks(),
    # xlabelvisible = false,
    # ygridvisible=false,
    # xgridvisible=false,
    # xticklabelsvisible= false,
    # xticksvisible= false,
    # # xminorticksvisible=false
    # ylabelvisible=false,
    # )

    # SqFTPl(x,y) = SqFieldTheory2(x,y,0/280,0.05/280)
    # fittingCoefs2 = fittingCoefs .* (10,0.5,0.5)
    SqFTPl = getFittedSq(fittingCoefs)
    SqFT_0 = SqFieldTheory

    SqFT0 = [SqFT_0(x,y) for x in kx, y in ky]
    SqFT = [SqFTPl(x,y) for x in kx, y in ky]
    normalizationFT0 = maximum(Sq)/maximum(SqFT0)
    SqFT0 .*= normalizationFT0
    # SqFT0 ./= maximum(SqFT0)
    # SqFT ./= maximum(SqFT)

    hmFT0 = heatmap!(axFT,kx,ky,SqFT0 ,colorrange = extrema(Sq),colormap = :viridis)

    hmFT = heatmap!(axFTB,kx,ky,SqFT ,colorrange = extrema(Sq),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))

    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)

    # kxB = kyB = LinRange(-pi,pi,200)
    kxB = kyB = trueMomenta(-pi,pi,size(BBMat,1)-1)
    BBPl = heatmap!(axBB,kxB,kyB,BBMat;
    colorrange = (min(0,minimum(BBMat)),maximum(BBMat)),colormap = Makie.cgrad(:thermal,rev=false)) # GnBu, dense, thermal, magma,plasma


    # path(t,angle) = t .*(cos(angle),sin(angle)) #.+ (pi,pi)
    
    
    # tRange = LinRange(0,.2pi,500)
    # angle = deg2rad(45)
    angle = pi/4
    # p1 = [path(t,angle) for t in tRange]

    
    xygrid = [(x,y) for x in kx, y in ky]

    tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
    # filter!(x -> x[1] in kx && x[2] in ky,p1)
    lines!(axFT,Point.(p1);color=:blue,linewidth = 1,linestyle = :dash)
    lines!(axFTB,Point.(p1);color=:red,linewidth = 1,linestyle = :solid)
    scatter!(axMC,Point.(xygrid[p1_discrete]);color=:black,markersize = 3)
    
    text!(axFT,Point(0,0),text="Γ",color = :white,align = (:center,:center))
    text!(axFT,Point(pi,0),text="X",color = :white,align = (:center,:center))
    text!(axFT,Point(0,pi),text="X'",color = :white,align = (:center,:center))
    # lines!(axFT,Point2f.([(0,0),(pi,0),(pi,pi),(0,pi),(0,0)]),color=:white,linewidth = 2,linestyle = :dash)
    # lines!(axFT,Point2f.([(-pi/2,pi/2),(pi/2,3pi/2),(3pi/2,pi/2),(pi/2,-pi/2),(-pi/2,pi/2)]),color=:white,linewidth = 2,linestyle = :dash)

    # lines!(axMC,Point.(p1);color=:black,linewidth = 2)

    colors = (:black,:green,:darkblue)
    markersizes = (4,4,4)
    # for (color,p) in zip(colors,(500,200))

    # axislegend(axpath1,merge=true,position=:lt)
    SqFTp = [SqFT_0(x,y) for (x,y) in (xygrid[p1_discrete])]
    # SqFTp = SqFT_0.(p1)
    SqFTp .*= normalizationFT0

    # SqFTp ./= maximum(SqFTp)
    lines!(axpath1,tRange,SqFTp;color=:blue,linewidth = 1.4,linestyle = :dash)    
    # lines!(axpath1,eachindex(p1).-1,SqFTp;color=:blue,linewidth = 1.4,linestyle = :dash)    

    SqFTp2 = [SqFTPl(x,y) for (x,y) in (xygrid[p1_discrete])]
    lines!(axpath1,tRange,SqFTp2;color=:red,linewidth = 1.6,linestyle = :solid)
    # SqFTp2 = SqFTPl.(p1)
    # lines!(axpath1,eachindex(p1).-1,SqFTp2;color=:red,linewidth = 1.6,linestyle = :solid)

    # SqFTp2 ./= maximum(SqFTp2)
    for (color,p,markersize) in zip(colors,(500,),markersizes)
        SqsGFMC_ps = eachslice(resOld["SqsGFMC"][string(p)],dims=3)./4
        SqsGFMC_p = mean(SqsGFMC_ps)
        # SqsGFMC = eachslice(res["30"]["SqsGFMC"][string(p)],dims=3)./4
        SqsGFMC_p_err = std(SqsGFMC_ps)
        
        Sf = SW.getSqCont(real(SqsGFMC_p))
        Sferr = SW.getSqCont(real(SqsGFMC_p_err))

        Sqp = [Sf(x,y) for (x,y) in (xygrid[p1_discrete])]
    
        Sqerr = [Sferr(x,y) for (x,y) in (xygrid[p1_discrete])]

        scatter!(axpath1,tRange,Sqp;color,markersize,label = L"\mathcal{P} = %$p")
        errorbars!(axpath1,tRange,Sqp,Sqerr;color,whiskerwidth = 4)
        # xlims!(axpath1,-0.1,0.85pi)
    end
    # axislegend(axpath1,merge=true,position=:lt)
    # lines!(axpath1,tRange²,tRange.^4*0.04/(0.03pi);color=:blue,linewidth = 2,linestyle = :solid)
    
    Colorbar(fig_top[1,4],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig_bottom_right[1,2],BBPl,height = Relative(0.8),ticks = SimpleTicks([0,0.5,1.0]),
    # label = L"\mathcal{B}^2(\textbf{q})",
    #vertical=false,#flip_vertical_label=false,flipaxis = false,
    )
    Label(fig_top[1,1, TopLeft()],L"a)$$",padding =(-45,0,-30,0))
    Label(fig_top[1,2, TopLeft()],L"b)$$",padding =(-20,0,-30,0))
    Label(fig_top[1,3, TopLeft()],L"c)$$",padding =(-20,0,-30,0))
    Label(fig_bottom_left[1,1, TopLeft()],L"d)$$",padding =(-45,10,10,0))
    Label(fig_bottom_right[1,1, TopLeft()],L"e)$$",padding =(-45,10,10,0))

    colgap!(fig_bottom_right,1,Relative(0.1))
    rowgap!(fig.layout,1,Relative(-0.1))
    # (fig.layout,1,Relative(-0.1))
    colsize!(fig_bottom_right,2,Relative(-0.2))
    # Box(fig_bottom_right[1,1])
    save("../figs/PaperFigs/Sq_L40_comp.pdf",fig)
    fig
end
##
fileBCorrRK = "../Data/BBCorrRK/BBCorrRK.h5"
BBCorrRK = h5open(fileBCorrRK) do f
    read(f)
end
##
function getMuSweep(BBCorrRK,p)
    BValsVec = Float64[]
    BValsErr = Float64[]
    mus = Float64[]
    # BValsErr = zeros(length(BBCorrRK["BVals"]))
    for (k,v) in BBCorrRK["BVals"]
        push!(BValsVec, SW.mean(v["B"],dims=2)[p])
        push!(BValsErr, SW.std(v["B"],dims=2)[p])
        push!(mus, v["mu"])
    end
    return mus,BValsVec,BValsErr
end
mus, BValsVec,BValsErr = getMuSweep(BBCorrRK,50)
##
with_theme(theme_SimpleTicks()) do 
    fig = Figure(size = (350,300))
    ax = Axis(fig[1,1],xlabel = L"μ",ylabel = L"B",aspect = 1)
    scatter!(ax,mus,BValsVec)
    errorbars!(ax,mus,BValsVec,BValsErr,whiskerwidth = 8)
    fig
end