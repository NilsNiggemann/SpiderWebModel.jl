import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using HDF5, H5Zblosc
##
files = readdir("../../Data/GFMC_1/",join=true)
AllResults = map(f -> h5open(f) do f
    @time read(f)
end,files)
##
function getEns(files)
    nThermal = 1_000
    ens = Vector{Vector{Float64}}()
    for file in files
        try
            h5open(file) do f
                en = SW.getEnergies(f["TotalWeights"][nThermal:end],f["energies"][nThermal:end],1,550÷nBra)
                push!(ens,en)
            end
        catch
        end
    end
    return ens
end
ens = getEns(files)

##
en = mean(ens)
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0/L^2",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))

    err = sqrt.(var(ens))
    proj = nBra .*eachindex(en)

    scatter!(ax,proj,en/40^2,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,proj,en/40^2,err/40^2,whiskerwidth = 3.5,color = :black)
    axislegend(ax,merge=true)
    # ylims!(ax,-75.05,-71.9)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end
##
@views function getSq(res,p,nThermal)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,p)
    # Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)
    Sq = similar(Conf, ComplexF64)
    
    Si = similar(Conf, ComplexF64)
    plan = SW.LatticeFFTs.FFTW.plan_fft(Conf)

    function SqFunc(Conf)
        Si .= Conf
        SW.mul!(Sq, plan, Si)
        Sq .= abs2.(Sq)
    end
    SaveConfs = res.SaveConfigs[:,:,:,nThermal:end]
    reconfTable = res.reconfTable[:,nThermal:end]
    res = SW.getObs(Gnp,SaveConfs,reconfTable,SqFunc,p÷2)
    newRes = similar(res,size(res).+1)
    newRes[begin:end-1,begin:end-1] .= res

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    newRes ./NSites
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getSqs(files,p)
    nThermal = 1_000
    Sqs = Vector{Matrix{Float64}}(undef,length(files))
    Threads.@threads for i in eachindex(files,Sqs)
        # file = files[i]
        f = files[i]
        # h5open(file) do f
            TotalWeights = f["TotalWeights"]
            SaveConfigs = f["SaveConfigs"]
            reconfTable = f["reconfigurationTable"]
            nBra = f["nBra"]
            res = (;TotalWeights,SaveConfigs,reconfTable)
            Sq = getSq(res,p÷nBra,nThermal)
            Sqs[i] = Sq
        # end
    end
    return Sqs
end
projectionSteps = 250
SqsGFMC = getSqs(AllResults,projectionSteps)
##

function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end
SqFieldTheory(k) = SqFieldTheory(k[1],k[2])
with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    SqMat = mean(real(SqsGFMC)) ./4
    SqFunc = SW.getSqCont(SqMat)
    
    # errMat = sqrt.(var(real(SqsGFMC))) ./4 ./SqMat
    errMat = sqrt.(var(real(SqsGFMC))) ./4
    errFunc = SW.getSqCont(errMat)
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .+ 0.5pi
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = (650,400))
    
    ticks = PiTicks()
    ticks = PiTicks([pi,2pi,])
    # ticks = PiTicks([-pi,-pi/2,0])
    # ticks = ([pi/2,3pi/2,5pi/2],[L"\frac{π}{2}",L"\frac{3}{2}π",L"\frac{5}{2}π"])
    # ticks = ([-3pi/2,-pi/2,pi/2],[L"-3π/2",L"-π/2",L"π/2"])
    # ticks = ([1,-1],[L"\frac{1}{2}",L"2"])
    axFT = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,xticks = ticks,yticks = ticks)
    axMC = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)
    axerr = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false,xticks = ticks,yticks = ticks)


    # axDiff = Axis(fig[1,4],xlabel = L"k_x",ylabel = L"k_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT ./maximum(SqFT),colormap = :viridis)
    # hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))
    hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(Sq))

    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)
    
    Colorbar(fig[2,1],hmFT,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8), width = Relative(0.9),vertical=false,ticks = SimpleTicks(),flipaxis=false)
    Colorbar(fig[2,2:3],hmMC,label = L"\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.9),ticks = SimpleTicks(),flipaxis=false)
    # Colorbar(fig[2,3],hmerr,label = L"\sigma\mathcal{S}^{zz}(\textbf{q})",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    rowgap!(fig.layout,1,-1.1)    
    colgap!(fig.layout,2,0.2)
    # colgap!(fig.layout,2,0.05)
    rowsize!(fig.layout,2,Relative(0.05))
    save("Application/exactFig/Spin1/Sq_L40.pdf",fig)
    fig
end
##

with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    SqMat = mean(real(SqsGFMC)) ./4
    SqFunc = SW.getSqCont(SqMat)
    
    # errMat = sqrt.(var(real(SqsGFMC))) ./4 ./SqMat
    errMat = sqrt.(var(real(SqsGFMC))) ./4
    errFunc = SW.getSqCont(errMat)
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .- 1.5pi
    kxSmall = kySmall = filter(x->abs(x) < pi/2,kx)

    err = [errFunc(x,y) for x in kxSmall, y in kySmall]

    fig = Figure(fontsize = 20,size = (600,400))
    
    ticks = ([-3pi/2,-pi/2,pi/2],[L"-3π/2",L"-π/2",L"π/2"])
    smallticks = ([-pi/3,0,pi/3],[L"-π/3",L"0",L"π/3"])

    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC $p=%$projectionSteps$",aspect = 1,xticks = ticks,yticks = ticks)
    
    
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std. error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=true,xticks = smallticks,yticks = smallticks)

    hmMC = halfhalfheatmap!(axMC,kx,ky,SqFieldTheory,SqFunc,x->-x-pi,normalize = true)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
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
    SqMat = mean(real(SqsGFMC)) ./4
    SqFunc = SW.getSqCont(SqMat)

    errMat = sqrt.(var(real(SqsGFMC))) ./4
    errFunc = SW.getSqCont(errMat)
    kx = ky = 2pi .* LinRange(0,1,size(SqMat,1)) .+ 0.5pi
    
    Sq = [SqFunc(x,y) for x in kx, y in ky]
    err = [errFunc(x,y) for x in kx, y in ky]

    fig = Figure(fontsize = 20,size = (650,450))
    
    ticks = PiTicks([pi,2pi,])

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
    save("Application/exactFig/Spin1/Sq_L40.pdf",fig)
    fig
end
##
