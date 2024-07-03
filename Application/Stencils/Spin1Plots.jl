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


##
with_theme(theme_PiTicks()) do
    Sq = SW.getSqCont(SqRK)
    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqRK,1)-1)
    heatmap(kx,ky,Sq.(Iterators.product(kx,ky)))
end
##
res = h5open("../Data/Spin1GFMC_Eval_periodic.h5") do f
    read(f)
end

projectionSteps = 500
SqsGFMC = eachslice(res["30"]["SqsGFMC"][string(projectionSteps)],dims=3)

SqMat = mean(real(SqsGFMC)) ./4
errMat = sqrt.(var(real(SqsGFMC))) ./4
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
    # save("Application/exactFig/GFMCEnergy.png",fig)
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
# function SqFieldTheory2(x,y,v,w)
#     a = (cos(x) - cos(y))^2 +(2sin(x)sin(y))^2
#     b = -((-4 + 4cos(x)*cos(y)+cos(2x)*(1-2cos(2y)) +cos(2y))*(v+8w-2w*(4cos(x)*cos(y)+cos(2x)*(1-2cos(y)) + cos(2y))))
#     c = (cos(x) - cos(y) + 2sin(x)sin(y))^2 
#     return a*√2 * sqrt(b)*c
# end
# function SqFieldTheory2(x,y,v,w)
#     kx =x 
#     ky = y
#     (sqrt(2)*sqrt(-((-4 + 4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky))*(v + 8*w - 2*w*(4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky)))))*(cos(kx) - cos(ky) + 2*sin(kx)*sin(ky))^2)/((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2+1e-2)
# end
# function SqFieldTheory2(kx,ky,v,w)
#     kx = kx -pi
#     ky = ky -pi
#     (sqrt(2)*sqrt(4 - 4*cos(kx)*cos(ky) - cos(2*ky) + cos(2*kx)*(-1 + 2*cos(2*ky)))*(-cos(kx) + cos(ky) + 2*sin(kx)*sin(ky))^2)/(sqrt(v + 8*w + 2*w*(-4*cos(kx)*cos(ky) - cos(2*ky) + cos(2*kx)*(-1 + 2*cos(2*ky))))*((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2) +1e-20)
# end

function SqFieldTheory2(kx,ky,k,b2,)
    (sqrt(2)*sqrt(-((-4 + 4*cos(kx)*cos(ky) + cos(2*kx)*(1 - 2*cos(2*ky)) + cos(2*ky))*(40*b2 + k + 8*b2*(cos(2*kx) + 2*cos(kx)*(-4*cos(ky) + cos(kx)*cos(2*ky))))))*   (cos(kx) - cos(ky) + 2*sin(kx)*sin(ky))^2)/((cos(kx) - cos(ky))^2 + 4*sin(kx)^2*sin(ky)^2+1e-30)
end

SqFieldTheory2(k,v,w) = SqFieldTheory2(k[1],k[2],v,w)

function optimizeCoeffs(SqMat)
    k = 2pi .* LinRange(0,1,size(SqMat,1))

    # SqFT = [SqFieldTheory2(x,y,0,0) for x in kx, y in ky]
    function loss(v,w)
        l = 0.
        v = abs(v)
        w = abs(w)
        for (i,kx) in enumerate(k), (j,ky) in enumerate(k)
            # l += abs2(SqMat[i,j] - SqFieldTheory2(kx,ky,v,w))#/(SqMat[i,j]+1e-5)
            # if (kx^2 +ky^2) < (0.8pi)^2
            l += abs2(SqMat[i,j] - SqFieldTheory2(kx,ky,v,w))#/(SqMat[i,j]+1e-5)
            # end
        end
        # for (i,kx) in enumerate(kx)
        #     l += abs2(SqMat[i,1] - SqFieldTheory2(kx,0,v,w))/(SqMat[i,1]+1e-5)
        # end
        return l
    end
    loss(v) = loss(v[1],v[2])

    # inner_optimizer = GradientDescent()
    # v0 = 0.08125
    # w0 = 0.025
    v0 = 1.
    w0 = 1.
    println(loss(v0,w0))
    x0 = [v0, w0]
    # (;v,w) = optimize(loss, [1e-10,1e-10],[Inf,Inf],Fminbox(inner_optimizer))
    res = optimize(loss, x0)
    @info res
    (v,w) = abs.(Optim.minimizer(res))
    @info "" loss(v0,w0) loss(v,w)
    return (;v,w)
end
fittingCoeffs::NamedTuple{(:v, :w), Tuple{Float64, Float64}} = optimizeCoeffs(SqMat)
##
SqFieldTheory(k) = SqFieldTheory2(k,fittingCoeffs...)
SqFieldTheory(kx,ky) = SqFieldTheory2(kx,ky,fittingCoeffs...)
##
# SqFieldTheory(k) = SqFieldTheory2(k,1/850,0.2/850)
# SqFieldTheory(kx,ky) = SqFieldTheory2(kx,ky,1/850,0.2/850)



with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 20,size = 2 .*(450,155))
    # fig = Figure()
    
    ticks = PiTicks([0,pi,])
    Ls = sort(collect(keys(res)))

    minval = Inf
    maxval = -Inf
    projectionSteps = 500

    for (L,resL) in res
        SqsGFMC = eachslice(resL["SqsGFMC"][string(projectionSteps)],dims=3)

        SqMat = mean(real(SqsGFMC)) ./4

        minval = min(minval,minimum(SqMat))
        maxval = max(maxval,maximum(SqMat))
    end

    for (i,Lstr) in enumerate(Ls)
        resL = res[Lstr]
        L = parse(Int,Lstr)
        SqsGFMC = eachslice(resL["SqsGFMC"][string(projectionSteps)],dims=3)

        SqMat = mean(real(SqsGFMC)) ./4
        SqMat ./= maxval

        # errMat = sqrt.(var(real(SqsGFMC))) ./4
        SqFunc = SW.getSqCont(SqMat)

        kx = ky = trueMomenta(-0.5pi,1.5pi,L)
        Sq = [SqFunc(x,y) for x in kx, y in ky]

        axFT = Axis(fig[1,i],xlabel = L"q_x",ylabel = L"q_y",title = L"$L=%$L$",aspect = 1,xticks = ticks,yticks = ticks,xlabelvisible=true,yticklabelsvisible=i==1,ylabelvisible = i==1,xticklabelsvisible=true)

        hmMC = halfhalfheatmap!(axFT,kx,ky,SqFieldTheory,SqFunc,x->-x+pi,normalize = true)

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
    save("Application/exactFig/Sq_size_scaling.pdf",fig)
    fig
end
##
with_theme(theme_SimpleTicks()) do 
    kx = trueMomenta(-0.5pi,1.5pi,30)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"\mathcal{S}(\mathbf{q})",xticks = PiTicks())
    
    colors = (:lightblue,:darkblue,:black)
    for (color,p) in zip(colors,(100,200,500))
        SqsGFMC = eachslice(res["30"]["SqsGFMC"][string(p)],dims=3)

        SqMat = mean(real(SqsGFMC)) ./4
        errMat = sqrt.(var(real(SqsGFMC))) ./4
        SqFunc = SW.getSqCont(SqMat)

        errFunc = SW.getSqCont(errMat)

        Sqerr = [errFunc(x,0) for x in kx]
        SqVec = [SqFunc(x,0) for x in kx]
        # SqVec = SqMat[1,:]
        # Sqerr = errMat[1,:]
        SqVec ./= maximum(SqVec)
        Sqerr ./= maximum(SqVec)

        scatterlines!(kx,SqVec;color)
        errorbars!(kx,SqVec,Sqerr;color,whiskerwidth = 8)
    end
    SqFT = [SqFieldTheory(x,0) for x in kx]
    SqFT ./= maximum(SqFT)

    lines!(kx,SqFT,color = :red)
    fig
end
##

with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = trueMomenta(-0.5pi,1.5pi,30)
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
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]

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
    # save("Application/exactFig/Spin1/Sq_L40.pdf",fig)
    fig
end
##
with_theme(theme_PiTicks()) do 

    kx = ky = trueMomenta(-0.5pi,1.5pi,size(SqMat,1)-1)
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
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
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]

    hmFT = heatmap!(axFT,kx,ky,SqFT ,colorrange = extrema(Sq),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    path(t,angle) = t .*(cos(angle),sin(angle)) #.+ (pi,pi)
    
    
    tRange = LinRange(0,.2pi,500)
    # angle = deg2rad(45)
    angle = pi/4
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
    SqFTp = [SqFieldTheory(x,y) for (x,y) in (xygrid[p1_discrete])]

    lines!(axpath2,tRange²,SqFTp;color=:red,linewidth = 2,linestyle = :dash)
    # lines!(axpath1,tRange²,tRange.^4*0.04/(0.03pi);color=:blue,linewidth = 2,linestyle = :solid)
    
    Colorbar(fig[1,3],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(1),ticks = SimpleTicks())
    Label(fig[1,1, TopLeft()],L"a)$$",padding =(-45,0,0,0))
    Label(fig[1,2, TopLeft()],L"b)$$",padding =(-20,0,0,0))
    Label(fig[2,1, TopLeft()],L"c)$$",padding =(-45,10,10,0))
    # rowgap!(fig.layout,1,Relative(-0.))
    rowgap!(fig.layout,1,Relative(-0.1))
    rowsize!(fig.layout,2,Relative(0.5))

    save("../exactFig/Spin1/Sq_L40_comp.pdf",fig)
    fig
end

##


with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = 2pi .* LinRange(0,1,500) .- 0.5pi
    
    fig = Figure(fontsize = 20,size = (450,400))
    
    ticks = PiTicks([0,pi,])
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks)

    SqFT = [SqFieldTheory2(x,y,5,1) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    
    Colorbar(fig[1,2],hmFT,label = L"\mathcal{S}(\textbf{q})",height = Relative(0.95),vertical=true,ticks = SimpleTicks(),flipaxis=true)
    # save("../Application/exactFig/Sq_fieldTheory.png",fig)
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
    save("Application/exactFig/omegaq_fieldTheory.png",fig)
    fig
end
##
BBCorrsFile = "../Data/BBCorr2/BBCorr_L=20.h5"
BBCorrs = h5read(BBCorrsFile,"BBCorrelator")
L = h5read(BBCorrsFile,"L")
nBra = only(unique!(h5read(BBCorrsFile,"nBra")))

BBCorrelator_p = mean(BBCorrs,dims=1)[1,:,:]
BBCorrelator_err_p = sqrt.(var(BBCorrs,dims=1))[1,:,:]
refPlaq = SW.getCentralPlaquette(S)
##
S = SW.stencilConfig(zeros(L,L),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
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
    save("Application/exactFig/Spin1/RealSpaceCorr_L20.pdf",fig)
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
    ri = [T * SW.SVector(r .- refPlaq) for r in AllPlaqs]
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
    save("Application/exactFig/Spin1/BCorrFT_L20_comp.pdf",fig)
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

with_theme(theme_PiTicks()) do
    T = SW.SA[
        1 1;
        -1 1
    ]/2
    ri = [T * SW.SVector(r .- refPlaq) for r in AllPlaqs]
    # rPlaq = [mapToPlaquetteBasis(r) for r in ri]
    capfilter(x) = min(abs(x),30)


    k = LinRange(-pi,pi,200)

    FTgen(BBCorrelator) = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    hm = nothing
    fig = Figure(size = 1. .*(400,320))
    p = 5
    FTs = FTgen.(eachslice(BBCorrs[:,:,p],dims=1))
    FT = mean(FTs)
    FTerr = sqrt.(var(FTs))#

    ax = Axis(fig[1,1];aspect = 1,ylabel = L"q_y",xlabel = L"q_x")
    # ax3 = Axis(fig[1,3];aspect = 1)
    
    hm = heatmap!(ax,k,k,FT)

    Colorbar(fig[1,2],hm,label = L"\mathcal{B}^2(\mathbf{q})",ticks = SimpleTicks())
    save("Application/exactFig/Spin1/BCorrFT_L20.pdf",fig)
    fig
end
