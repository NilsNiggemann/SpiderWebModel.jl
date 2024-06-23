using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW
using Optim
cd(@__DIR__)
##
filesRK = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/DataRK/",join=true)
resRK = h5open(filesRK[end]) do f
    read(f)
end
##
res = h5open("../Data/Spin1GFMC_Eval_open.h5") do f
    read(f)
end
##
with_theme(theme_SimpleTicks()) do
    en = res["energies"]
    e0s,e0s_err = mean(en,dims=2)[:],sqrt.(var(en,dims=2)[:])
    # e0s = (e0s .- minimum(e0s) .+1e-1)
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0/L^2",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))

    proj = res["nBra"] .*eachindex(e0s)
    scatter!(ax,proj,e0s/40^2,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,proj,e0s/40^2,e0s_err/40^2,whiskerwidth = 3.5,color = :black)
    axislegend(ax,merge=true)
    # ylims!(ax,-75.05,-71.9)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end
##
projectionSteps = 1000
SqsGFMC = eachslice(res["SqsGFMC"][string(projectionSteps)],dims=3)

SqMat = mean(real(SqsGFMC)) ./4
errMat = sqrt.(var(real(SqsGFMC))) ./4
SqFunc = SW.getSqCont(SqMat)

errFunc = SW.getSqCont(errMat)

# function SqFieldTheory(x,y)
#     num = cos(x) - cos(y) +2sin(x)sin(y) 
#     denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
#     return num^2/(sqrt(denom)+1e-30)
# end
# SqFieldTheory(k) = SqFieldTheory(k[1],k[2])
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
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = (650,400))
    
    ticks = PiTicks([0,pi,])
    # ticks = PiTicks([-pi,-pi/2,0])
    # ticks = ([pi/2,3pi/2,5pi/2],[L"\frac{π}{2}",L"\frac{3}{2}π",L"\frac{5}{2}π"])
    # ticks = ([-3pi/2,-pi/2,pi/2],[L"-3π/2",L"-π/2",L"π/2"])
    # ticks = ([1,-1],[L"\frac{1}{2}",L"2"])
    axFT = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"q_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks)
    axMC = Axis(fig[1,2],xlabel = L"q_x",ylabel = L"q_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)
    axerr = Axis(fig[1,3],xlabel = L"q_x",ylabel = L"q_y",title = L"rel. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)


    # axDiff = Axis(fig[1,4],xlabel = L"q_x",ylabel = L"q_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    hmerr = heatmap!(axerr,kx,ky,err ./Sq,colormap = :viridis)

    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)
    
    Colorbar(fig[2,1],hmFT,label = L"\mathcal{S}(\textbf{q})",height = Relative(0.8), width = Relative(0.9),vertical=false,ticks = SimpleTicks(),flipaxis=false)
    Colorbar(fig[2,2],hmMC,label = L"\mathcal{S}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.9),ticks = SimpleTicks(),flipaxis=false)
    
    Colorbar(fig[2,3],hmerr,label = L"\sigma\mathcal{S}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks(),flipaxis=false)
    rowgap!(fig.layout,1,-1.1)    
    colgap!(fig.layout,2,0.2)
    # colgap!(fig.layout,2,0.05)
    rowsize!(fig.layout,2,Relative(0.05))
    # save("Application/exactFig/Spin1/Sq_L40.pdf",fig)
    fig
end
##
with_theme(theme_SimpleTicks()) do 
    kx = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"q_x",ylabel = L"\mathcal{S}(\mathbf{q})",xticks = PiTicks())
    
    colors = (:lightblue,:darkblue,:black)
    for (color,p) in zip(colors,(500,750,1000))
        SqsGFMC = eachslice(res["SqsGFMC"][string(p)],dims=3)

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
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
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

    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
    
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

    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
    
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
    for (color,p) in zip(colors,(1000,500))
        SqsGFMC = eachslice(res["SqsGFMC"][string(p)],dims=3)./4

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
    lines!(axpath1,tRange²,tRange.^4*0.04/(0.03pi);color=:blue,linewidth = 2,linestyle = :solid)
    
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