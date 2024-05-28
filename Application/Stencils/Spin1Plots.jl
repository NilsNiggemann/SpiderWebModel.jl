using CairoMakie, MakieHelpers,Statistics, HDF5
import SpiderWebModel as SW

##
filesRK = readdir("/p/scratch/pmfrg/niggemann1/Spiderweb/DataRK/",join=true)
resRK = h5open(filesRK[1]) do f
    read(f)
end
##
res = h5open("../Data/Spin1GFMC_Eval.h5") do f
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
projectionSteps = 250
SqsGFMC = eachslice(res["SqsGFMC"][string(projectionSteps)],dims=3)

SqMat = mean(real(SqsGFMC)) ./4
errMat = sqrt.(var(real(SqsGFMC))) ./4
SqFunc = SW.getSqCont(SqMat)

errFunc = SW.getSqCont(errMat)

function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end
SqFieldTheory(k) = SqFieldTheory(k[1],k[2])
##


with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 0.5pi
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = (650,400))
    
    ticks = PiTicks()
    ticks = PiTicks([0,pi,])
    # ticks = PiTicks([-pi,-pi/2,0])
    # ticks = ([pi/2,3pi/2,5pi/2],[L"\frac{π}{2}",L"\frac{3}{2}π",L"\frac{5}{2}π"])
    # ticks = ([-3pi/2,-pi/2,pi/2],[L"-3π/2",L"-π/2",L"π/2"])
    # ticks = ([1,-1],[L"\frac{1}{2}",L"2"])
    axFT = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks)
    axMC = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)
    axerr = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"rel. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)


    # axDiff = Axis(fig[1,4],xlabel = L"k_x",ylabel = L"k_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    hmerr = heatmap!(axerr,kx,ky,err ./Sq,colormap = :viridis)

    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)
    
    Colorbar(fig[2,1],hmFT,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.9),vertical=false,ticks = SimpleTicks(),flipaxis=false)
    Colorbar(fig[2,2],hmMC,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.9),ticks = SimpleTicks(),flipaxis=false)
    
    Colorbar(fig[2,3],hmerr,label = L"\sigma\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks(),flipaxis=false)
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
    SqFT = [SqFieldTheory(x,0) for x in kx]
    Sqerr = [errFunc(x,0) for x in kx]
    SqVec = [SqFunc(x,0) for x in kx]
    # SqVec = SqMat[1,:]
    # Sqerr = errMat[1,:]
    SqFT ./= maximum(SqFT)
    SqVec ./= maximum(SqVec)
    Sqerr ./= maximum(SqVec)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"\mathcal{S}(\mathbf{q})",xticks = PiTicks())
    scatterlines!(kx,SqVec,color = :black)
    errorbars!(kx,SqVec,Sqerr,color = :black,whiskerwidth = 8)
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

    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,xticks = ticks,yticks = ticks)
    
    
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=true,xticks = smallticks,yticks = smallticks)

    hmMC = halfhalfheatmap!(axMC,kx,ky,SqFieldTheory,SqFunc,x->-x-pi,normalize = true)
    # SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    # heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(!isnan,Sq))
    hmerr = heatmap!(axerr,kxSmall,kySmall,err,colormap = :viridis)
    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)

    Colorbar(fig[2,1],hmMC,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,2],hmerr,label = L"\sigma\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    
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

    axFT = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    axMC = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    # axerr = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    
    axpath = Axis(fig[2,1:3],xlabel = L"t",ylabel = L"\tilde{\mathcal{S}}^{zz}(\mathbf{q})",yticks = SimpleTicks(),yminorticksvisible = true,xminorticksvisible = true)

    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]

    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))

    tRange = LinRange(0,2pi,100)
    path(r,t) = (r*cos(t)+pi,r*sin(t)+pi)

    p1 = [path(0.65,t) for t in tRange]
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
    scatterlines!(axpath,tRange,SqFTp,color = :red,linwidth = 1,linestyle = :solid)
    scatter!(axpath,tRange,Sqp,color = :black,markersize = 8)
    errorbars!(axpath,tRange,Sqp,Sqerr,color = :black,whiskerwidth = 8)
    Colorbar(fig[1,3],hmMC,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(1),ticks = SimpleTicks())
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

    fig = Figure(fontsize = 20,size = (480,450))
    
    ticks = PiTicks([0,pi,])

    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"QMC$$",aspect = 1,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    axFT = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks,yminorticksvisible = true,xminorticksvisible = true)
    
    axpath1 = Axis(fig[2,1:3],xlabel = L"t²",ylabel = L"\tilde{\mathcal{S}}^{zz}(\mathbf{q})",yticks = SimpleTicks(),xticks = SimpleTicks(),yminorticksvisible = true,xminorticksvisible = true,xminorgridvisible =true,yminorgridvisible=true,yaxisposition=:left)
    
    axpath2 = Axis(fig[2,1:3],xlabel = L"t²",yticks = SimpleTicks(),xticks = SimpleTicks(),
    yaxisposition=:right,yticklabelcolor=:red,
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

    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    path(t,angle) = t .*(cos(angle),sin(angle)) #.+ (pi,pi)
    
    
    colors = (:black, :red)
    for (color,angle,axpath) in zip(colors, (pi/3,0.6pi),(axpath1,axpath2))
        
        tRange = LinRange(0.,0.5pi,500)
        p1 = [path(t,angle) for t in tRange]
        xygrid = [(x+0.5pi,y+0.5pi) for x in kx, y in ky]
        tRange,p1_discrete = rasterCurve(p1,xygrid,tRange)
        tRange² = (tRange).^2
        # filter!(x -> x[1] in kx && x[2] in ky,p1)
        lines!(axFT,Point.(xygrid[p1_discrete]);color,linewidth = 2)
        lines!(axMC,Point.(xygrid[p1_discrete]);color,linewidth = 2)
        Sqp_all = [S[p1_discrete] for S in real(SqsGFMC)]
        for S in Sqp_all
            S ./= maximum(S)
        end
        
        Sqp = mean(Sqp_all)
        @info "" first(p1_discrete) first(Sqp)
        
        # Sqp ./= maximum(Sq)
        # SqFT_coarse = SW.getSqCont(SqFT)
        Sqerr = sqrt.(var(Sqp_all))
        # SqFTp = SqFieldTheory.(p1)
        # SqFTp = SqFT[p1_discrete]
        SqFTp = SqFieldTheory.(xygrid[p1_discrete])
        SqFTp ./= maximum(SqFTp)

        lines!(axpath,tRange²,SqFTp;color,linwidth = 1,linestyle = :dash)
        scatter!(axpath,tRange²,Sqp;color,markersize = 8)
        errorbars!(axpath,tRange²,Sqp,Sqerr;color,whiskerwidth = 8)
        xlims!(axpath,-0.1,0.85pi)
    end
    Colorbar(fig[1,3],hmMC,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(1),ticks = SimpleTicks())
    rowgap!(fig.layout,1,Relative(-0.))
    rowsize!(fig.layout,2,Relative(0.45))
    # save("Application/exactFig/Spin1/Sq_L40.pdf",fig)
    fig
end
##
