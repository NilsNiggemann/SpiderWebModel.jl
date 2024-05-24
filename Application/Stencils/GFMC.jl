import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
##
S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# S = SW.stencilConfig(SW.constructConfigPath(15,15,SW.ALLGS_S12),1/2)
HStair = SW.generateHilbertSpace(SW.SpinConfig(S))
HTest = Array(HStair.H)
for i in 1:size(HTest,1)
    SW.normalize!(@view(HTest[i,:]),1)
end
ExSol = SW.SolveHKrylov(HTest)
E0 = ExSol.values[1]

getRKWavefunction(ψ) = SW.normalize!(one.(ψ))

v0 = ExSol.vectors[1]
# v0 = getRKWavefunction(ExSol.vectors[1])
HConfs = SW.spinConfig.(HStair.AllStates,Ref(SW.
SpinConfig(S)),Ref(HStair.plaqMapping))
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf)
        ind = AllSTDict[Conf]
        return v0[ind] + 0.01
    end
end

function makeBlocks(arr;blocksize = nothing, numBlocks = 32)
    if blocksize === nothing
        blocksize = length(arr) ÷ numBlocks
    end
    blockedArr = collect(SW.splitIntoBins(arr,blocksize))
    
    minsize,maxsize = extrema(length.(blockedArr))
    
    if maxsize!=minsize
        pop!(blockedArr)
    end
    return blockedArr
end
function plotEnergies(results,nBra,E0;Emin=E0-1e-2,Emax=E0+2e-2)
    ens = [SW.getEnergies(res.TotalWeights,res.energies,1,250÷nBra) for res in results]

    en = mean(ens)
    # errs = getErrBlocking(results[1].energies[nThermal:end],results[1].TotalWeights[nThermal:end],2*10^4,20,E0) ./ ( (length(results[1].energies)-nThermal) ÷ 2*10^4)
    ##
    with_theme(theme_SimpleTicks()) do
        fig = Figure(fontsize = 22)
        ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
        # ens = getfield.(obs,:E0)
        en = mean(ens)
        # M = length(results[1].energies)
        # Mk = M ÷ length(ens)
        # println(Mk)
        err = sqrt.(var(ens))
        # err = 0.004 .* ones(length(en))
        proj = nBra .*eachindex(en)
        scatter!(ax,proj,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
        errorbars!(ax,proj,en,err,whiskerwidth = 3.5,color = :black)
        !isnan(E0)&& hlines!([E0],color = :red,label = L"exact $$")
        axislegend(ax,merge=true)
        xlims!(ax,0.5,last(proj))
        !isnan(E0)&& ylims!(ax,Emin,Emax)
        # save("Application/exactFig/GFMCEnergy.png",fig)
        fig
    end
end
magEx = SW.getMagnetization(HConfs, v0)
##
@time SqEx = SW.getStructureFac(HConfs,v0)
# @time SqEx = SW.getEqualWeightStructureFac(HConfs)
# @time SqEx = SW.getStructureFac(HConfs,SW.normalize!(ones(length(v0))))
#___________ManyWalkers_______________________
##
# S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# ψG = SW.PlaquetteNumberGuidingFunction(0.197)
ψG = SW.PlaquetteNumberGuidingFunction(0.0)

nThermal = 10_000
SW.Random.seed!(1234)
# results = [SW.startManyWalkerGFMC(S,2,55_000,3,nThermal,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:35]
nBra = 1
##
# outfile = "Data/temp/S12/"
# rm(outfile;recursive=true,force=true)
# mkpath(outfile)
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,16,50000÷nBra,nBra,ψG,1;equilibration_steps=nThermal,outfile =string(outfile,i,".h5")) for i in 1:8])
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,1,30000÷nBra,nBra,ψG,10;equilibration_steps=nThermal) for i in 1:32])
# @time SW.startManyWalkerGFMC(S,16,nThermal+500_00÷nBra,nBra,ψG,1;outfile = string(outfile,1,".h5"))
##

# nBra = 1
# @time results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1200_000,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:6*4])
# ens = [SW.getEnergies(res.TotalWeights[nThermal:end],res.energies[nThermal:end],1,150÷nBra) for res in results]

plotEnergies(results,nBra,E0)
# plotEnergies(results,nBra,E0)
##
#___________Observables_______________________
##
# getAvgMag(Conf) = sum(Conf) ./ (2*length(Conf))
function getMag(results,pmax,I)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,pmax) for res in results]
    getMag(Conf) = Conf[I] /2

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs,res.reconfigurationTable,getMag,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getMagRK(results,p,I)
    # I = CartesianIndex(I)
    TotWeights = [SW.mean(res.TotalWeights) for res in results]
    weights = [abs2.(res.TotalWeights ./ W) for (res,W) in zip(results,TotWeights)]
    # SW.normalize!.(weights,1)
    i,j = Tuple(I)
    @views getObs(ws,res) = SW.mean(w * mag/2 for (w,mag) in zip(ws,res.SaveConfigs[i,j,1,:]))

    obs = fetch.([Threads.@spawn getObs(ws,res) for (ws,res) in zip(weights,results)])
    return [obs for i in 1:p]
end

function getSiSj(results,pmax,I,J)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,pmax) for res in results]

    getCorr(Conf) = (Conf[I]*Conf[J]) *0.25

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs,res.reconfigurationTable,getCorr,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
##
mags1 = getMagRK(results,100÷nBra,CartesianIndex(3,3)) 
mags2 = getMagRK(results,100÷nBra,CartesianIndex(4,1)) 
mags3 = getMagRK(results,100÷nBra,CartesianIndex(2,3)) 
##
IJ_SS = (CartesianIndex(3,4),CartesianIndex(4,2))
SiSj = getSiSj(results,100÷nBra,IJ_SS[1],IJ_SS[2]) 
SiSjex = SW.getSij(HConfs,v0,IJ_SS[1],IJ_SS[2])
##

function errorBarLegend(size = 0.5;linekwargs = (;),markerkwargs = (;))
    center = 0.5
    ymin = center - size/2
    ymax = center + size/2

    errorbar = [Point2f(center, ymin), Point2f(center, ymax)]
    [
        LineElement(linepoints = errorbar),
        MarkerElement(points = errorbar, marker = :hline, markersize = 10;linekwargs...),
        MarkerElement(points = [Point2f(center, center),], marker = '●', markersize = 7;markerkwargs...)
    ]
end

with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    proj = nBra .*eachindex(mags1)
    
    Cols = (:blue,:green,:orange)
    Is = ((3,3),(4,1),(2,3))
    for (i,mags) in enumerate((mags1,mags2,mags3))
        mExact = magEx[Is[i]...]

        sgn = sign(mExact)
        m =  sgn .* mean.(mags)
        scatter!(ax,proj,m,label = L"GFMC$$",color = Cols[i], marker = '●',markersize = 5)
        err =sqrt.(var.(mags))
        errorbars!(ax,proj,m,err,whiskerwidth = 3.5,color = Cols[i])
        
        hlines!([sgn * mExact],color = Cols[i],label = L"exact $$")
        sgstring = sgn > 0 ? "" : "-"
        text!(ax,(proj[end÷2],mExact-0.015),color = Cols[i],text = L"%$sgstring \langle S^z_{%$(Is[i])}\rangle")
    end
    chi = mean.(SiSj)

    chiScale = -3
    # scatter!(ax,proj,chiScale .*chi,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    # err =sqrt.(var.(SiSj))
    # errorbars!(ax,proj,chiScale .*chi,chiScale .*err,whiskerwidth = 3.5,color =  :black)
    # hlines!([chiScale*SiSjex],color = :black,label = L"exact $$")
    i,j = Tuple.(IJ_SS)
    text!(ax,(proj[end÷2],chiScale*SiSjex-0.015),color = :black,text = L"%$chiScale\times \langle S^z_{%$i} S^z_{%$j}\rangle",align = (:center,:center))


    Legend(fig[1, 1], [[
        errorBarLegend(0.6)
    ],
    [LineElement(linepoints = [Point2f(0., 0.5), Point2f(1, 0.5)])]], [L"GFMC$$",L"exact $$"], tellheight = false, tellwidth = false,halign = :right, valign = :center,margin = (10,10,120,10))


    # axislegend(ax,merge=true,position = :rc,unique=true)
    # xlims!(ax,0.5,last(proj))
    # ylims!(ax,E0-1e-2,E0+3e-2)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end
##
@views function getFullMag(res,p)
    # I = CartesianIndex(I)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    buffer = zeros(size(res.SaveConfigs[:,:,begin,begin]))
    function m(Conf)
        buffer .= Conf
    end
    res = SW.getObs(Gnp,res.SaveConfigs,res.reconfigurationTable,m,p÷2)
end
magFull = fetch.([Threads.@spawn getFullMag(res,100÷nBra) for res in results])

##
with_theme(theme_SimpleTicks()) do 
    m = mean(magFull) ./2
    fig,ax,hm = heatmap(m,colormap = :greys,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = 0.8 .*(360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,magEx,colormap = :greys,colorrange = extrema(m))
    Colorbar(fig[1:2,2],hm,label = L"\langle S^z_i \rangle")
    fig
end
##
function getSq(res,p)
    @views Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    # @views Gnp = ones(length(res.TotalWeights[nThermal:end]),p)

    Conf = @view res.SaveConfigs[:,:,begin,begin]
    NSites = length(Conf)
    Sq = similar(Conf, ComplexF64)
    
    Si = similar(Conf, ComplexF64)
    plan = SW.LatticeFFTs.FFTW.plan_fft(Conf)

    function SqFunc(Conf)
        Si .= Conf
        SW.mul!(Sq, plan, Si)
        Sq .= abs2.(Sq)
    end

    @views res = SW.getObs(Gnp,res.SaveConfigs,res.reconfigurationTable,SqFunc,p÷2)
    newRes = similar(res,size(res).+1)
    newRes[begin:end-1,begin:end-1] .= res

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    newRes ./NSites
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
##
SqsGFMC = fetch.([Threads.@spawn getSq(res,100÷nBra) for res in results])
##
with_theme(theme_PiTicks()) do 
    Sq = mean(real(SqsGFMC)) ./4
    (;kx,ky) = Sq
    fig,ax,hm = heatmap(kx,ky,Sq,colormap = :viridis,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = (360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,kx,ky,real(SqEx.Sq),colormap = :viridis,colorrange = extrema(Sq))
    Colorbar(fig[1:2,2],hm,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle")
    fig
end
##
#___________Spin-1_______________________
SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(15,15),1)
nBra = 5
nThermal = 0
##
function variationalEnergy(res)
    Es = mean(res.energies for res in res)
    
    E = mean(Es)
    err = sqrt(var(Es))
    return E,err
end
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,ψG,1) for _ in 1:8])
function getVarEnergies(alphaRange;nBra,nSteps,nWalkers,nRuns = 6)
    Es = zeros(size(alphaRange))
    errs = zeros(size(alphaRange))

    for (i,α) in enumerate(alphaRange)
        ψG = SW.PlaquetteNumberGuidingFunction(α)
        @time results = [SW.startManyWalkerGFMC(S,nWalkers,nSteps÷nBra,nBra,ψG,1) for _ in 1:nRuns]
        E,err = variationalEnergy(results)
        println("α = $α, E = $E, err = $err")
        Es[i] = E
        errs[i] = err
    end
    return Es,errs
end
alpharange = LinRange(0.1,0.16,20)
E,err = getVarEnergies(alpharange;nBra = 2,nSteps = 12_000,nWalkers = 4,nRuns = 24)
##
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"\alpha",ylabel = L"E^\textrm{var}_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))
    scatter!(ax,alpharange,E,label = L"$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,alpharange,E,err,whiskerwidth = 3.5,color = :black)
    fig

    
end
##
SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(16,16),1;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
nBra = 10
nThermal = 10_000
ψG = SW.PlaquetteNumberGuidingFunction(0.17)
##
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,10,120_000÷nBra,nBra,ψG,1;equilibration_steps=nThermal) for _ in 1:6])

##
plotEnergies(results,nBra,NaN)

##
SqsGFMC = fetch.([Threads.@spawn getSq2(res,80÷nBra) for res in results[1:6]])
##
function SqFieldTheory(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end

with_theme(theme_PiTicks()) do 
    # Sq = sqrt.(var(real(SqsGFMC))) ./4
    Sq = mean(real(SqsGFMC)) ./4
    kx = ky = 2pi .* LinRange(0,1,size(Sq,1))
    fig = Figure(fontsize = 22,size = (800,400))
    axMC = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"k_y",title = L"GFMC$$",aspect = 1)
    axerr = Axis(fig[1,2],xlabel = L"k_x",ylabel = L"k_y",title = L"std error$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    axFT = Axis(fig[1,3],xlabel = L"k_x",ylabel = L"k_y",title = L"U(1) theory$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    # axDiff = Axis(fig[1,4],xlabel = L"k_x",ylabel = L"k_y",title = L"Difference$$",aspect = 1,ylabelvisible = false,yticklabelsvisible=false)

    # err = abs.(Sq .- SqFieldTheory.(kx,kx'))
    err = sqrt.(var(real(SqsGFMC))) ./4
    hmMC = heatmap!(axMC,kx,ky,Sq,colormap = :viridis)
    SqFT = [SqFieldTheory(x,y) for x in kx, y in ky]
    hmFT = heatmap!(axFT,kx,ky,SqFT,colormap = :viridis)
    # heatmap!(axerr,kx,ky,err,colormap = :viridis,colorrange = extrema(!isnan,Sq))
    hmerr = heatmap!(axerr,kx,ky,err,colormap = :viridis)
    # heatmap!(axDiff,kx,ky,(Sq ./maximum(Sq)) .- (SqFT ./maximum(SqFT)),colormap = :viridis)

    Colorbar(fig[2,1],hmMC,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,2],hmerr,label = L"\sigma(\langle \mathcal{S}^{zz}(\textbf{q})\rangle)",height = Relative(0.8),vertical=false,width = Relative(0.8),ticks = SimpleTicks())
    Colorbar(fig[2,3],hmFT,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle",height = Relative(0.8), width = Relative(0.8),vertical=false,ticks = SimpleTicks())

    rowsize!(fig.layout,2,Relative(0.1))
    fig
end
##
projection_orders = 150:-50:50
SqsGFMC_p = [ fetch.([Threads.@spawn getSq(res,p÷nBra) for res in results]) for p in projection_orders]
##
# with_theme(merge(theme_PiTicks(),theme_dark())) do 
with_theme(theme_PiTicks()) do 
    fig = Figure(fontsize = 18,size = (500,400))
    ax = Axis(fig[1,1],xlabel = L"k_x",ylabel = L"\mathcal{S}(\mathbf{q})",title = L"GFMC$$",yticks = SimpleTicks())
    kx = 2pi .* LinRange(0,1,size(S,1).+1)

    for (p,Sqs) in zip(projection_orders,SqsGFMC_p)
        Sq = mean(real(Sqs))[:,1] ./4
        Sqerr = sqrt.(var(real(Sqs)))[:,1] ./4
        Sq ./= maximum(Sq)
        Sqerr ./= maximum(Sq)
        lines!(ax,kx,Sq,color = p,colorrange = extrema(projection_orders),colormap = :viridis,linewidth = 3,label = L"$p=%$(p)$")
        errorbars!(ax,kx,Sq,Sqerr,whiskerwidth = 3.5)
    end

    kx = 2pi .* LinRange(0,1,300)
    SqFT = SqFieldTheory.(kx,0)
    SqFT ./= maximum(SqFT)
    lines!(ax,kx,SqFT,color = :black,linewidth = 2,linestyle = :dash,label = L"U(1)$$")
    axislegend(ax,merge=true)
    fig
end