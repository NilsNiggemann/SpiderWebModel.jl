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
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1]
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
function plotEnergies(results,nBra,nThermal,E0;Emin=E0-1e-2,Emax=E0+2e-2)
    ens = [SW.getEnergies(w,e,1,250÷nBra) for res in results[1:end] for (w,e) in zip(makeBlocks(res.TotalWeights[nThermal:end],numBlocks=1),makeBlocks(res.energies[nThermal:end],numBlocks=1)) ]

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
        hlines!([E0],color = :red,label = L"exact $$")
        axislegend(ax,merge=true)
        xlims!(ax,0.5,last(proj))
        ylims!(ax,Emin,Emax)
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
ψG = SW.PlaquetteNumberGuidingFunction(0.197)

nThermal = 6_000
SW.Random.seed!(1234)
# results = [SW.startManyWalkerGFMC(S,2,55_000,3,nThermal,SW.ConstructVaritationalFunc(0.197,S),0) for _ in 1:35]
nBra = 6
##
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,16,nThermal+300_000÷nBra,nBra,ψG,1) for _ in 1:8])
##
# nBra = 1
# @time results = fetch.([Threads.@spawn SW.startSingleWalkerGFMC(S,nThermal+1200_000,SW.ConstructVaritationalFunc(0.197,S),1) for _ in 1:6*4])
# ens = [SW.getEnergies(res.TotalWeights[nThermal:end],res.energies[nThermal:end],1,150÷nBra) for res in results]

plotEnergies(results,nBra,nThermal,E0)
##
#___________Observables_______________________
##
# getAvgMag(Conf) = sum(Conf) ./ (2*length(Conf))
function getMag(results,pmax,I,nThermal)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,pmax) for res in results]

    getMag(Conf) = Conf[I] /2

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs[:,:,:,nThermal:end],res.reconfTable[:,nThermal:end],getMag,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end

function getSiSj(results,pmax,I,J,nThermal)
    # I = CartesianIndex(I)
    Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,pmax) for res in results]

    getCorr(Conf) = (Conf[I]*Conf[J]) *0.25

    @views getObs(p) = [SW.getObs(Gnp[:,1:p],res.SaveConfigs[:,:,:,nThermal:end],res.reconfTable[:,nThermal:end],getCorr,p÷2) for (res,Gnp) in zip(results,Gnps)]
    obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
##
mags1 = getMag(results,100÷nBra,CartesianIndex(3,3),nThermal) 
mags2 = getMag(results,100÷nBra,CartesianIndex(4,1),nThermal) 
mags3 = getMag(results,100÷nBra,CartesianIndex(2,3),nThermal) 
##
IJ_SS = (CartesianIndex(3,4),CartesianIndex(4,2))
SiSj = getSiSj(results,100÷nBra,IJ_SS[1],IJ_SS[2],nThermal) 
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
    scatter!(ax,proj,chiScale .*chi,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    err =sqrt.(var.(SiSj))
    errorbars!(ax,proj,chiScale .*chi,chiScale .*err,whiskerwidth = 3.5,color =  :black)
    hlines!([chiScale*SiSjex],color = :black,label = L"exact $$")
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
@views function getFullMag(res,p,nThermal)
    # I = CartesianIndex(I)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,p)
    buffer = zeros(size(res.SaveConfigs[1][1]))
    function m(Conf)
        buffer .= Conf ./2
    end
    res = SW.getObs(Gnp,res.SaveConfigs[:,:,:,nThermal:end],res.reconfTable[:,nThermal:end],m,p÷2)
end
magFull = fetch.([Threads.@spawn getFullMag(res,100÷nBra,nThermal) for res in results])

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
function getSq(res,p,nThermal)
    @views Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,p)
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

    @views res = SW.getObs(Gnp,res.SaveConfigs[:,:,:,nThermal:end],res.reconfTable[:,nThermal:end],SqFunc,p÷2)
    newRes = similar(res,size(res).+1)
    newRes[begin:end-1,begin:end-1] .= res

    @views newRes[end,begin:end] .= newRes[begin,:]
    @views newRes[begin:end,end] .= newRes[:,begin]
    newRes ./NSites
    # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
end
# @views function getSq(res,p,nThermal)
#     # I = CartesianIndex(I)
#     Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights[nThermal:end],1,p)

#     Sq(x) = abs2.(SW.LatticeFFTs.fft(x))

#     res = SW.getObs(Gnp,res.SaveConfigs[:,:,:,nThermal:end],res.reconfTable[:,nThermal:end],Sq,p÷2)
#     newRes = similar(res,size(res).+1)
#     newRes[1:end-1,1:end-1] .= res
#     newRes[end,1:end-1] .= res[1,:]
#     newRes[1:end-1,end] .= res[:,1]
#     newRes./length(res)
#     # obs = fetch.([Threads.@spawn getObs(p) for p in 1:pmax])
# end
##
SqsGFMC = fetch.([Threads.@spawn getSq(res,100÷nBra,nThermal) for res in results])
##
with_theme(theme_PiTicks()) do 
    Sq = mean(real(SqsGFMC)) ./4
    (;kx,ky) = SqEx
    fig,ax,hm = heatmap(kx,ky,Sq,colormap = :viridis,axis=(;aspect=1,title = L"GFMC$$"),figure = (;size = (360,500)))
    ax2 = Axis(fig[2,1],aspect=1,title = L"exact $$")
    heatmap!(ax2,kx,ky,real(SqEx.Sq),colormap = :viridis,colorrange = extrema(Sq))
    Colorbar(fig[1:2,2],hm,label = L"\langle \mathcal{S}^{zz}(\textbf{q})\rangle")
    fig
end
##
#___________Spin-1_______________________
SW.Random.seed!(1234)
S = SW.stencilConfig(zeros(40,40),1)
nBra = 10
ψG = SW.PlaquetteNumberGuidingFunction(0.15) 
nThermal = 1
##
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,ψG,1) for _ in 1:8])
@time results = [SW.startManyWalkerGFMC(S,8,nThermal+3_000÷nBra,nBra,ψG,1) for _ in 1:8]
##
ens = [SW.getEnergies(w,e,1,200÷nBra) for res in results[1:end] for (w,e) in zip(makeBlocks(res.TotalWeights[nThermal:end],numBlocks=1),makeBlocks(res.energies[nThermal:end],numBlocks=1)) ]

en = mean(ens)
with_theme(theme_SimpleTicks()) do
    fig = Figure(fontsize = 22)
    ax = Axis(fig[1,1],xlabel = L"projection order $$",ylabel = L"E_0",xminorticksvisible=true,yminorticksvisible=true,xminorticks=IntervalsBetween(5),yminorticks = IntervalsBetween(5))

    err = sqrt.(var(ens))
    proj = nBra .*eachindex(en)

    scatter!(ax,proj,en,label = L"GFMC$$",color = :black, marker = '●',markersize = 5)
    errorbars!(ax,proj,en,err,whiskerwidth = 3.5,color = :black)
    axislegend(ax,merge=true)
    # ylims!(ax,-75.05,-71.9)
    # save("Application/exactFig/GFMCEnergy.png",fig)
    fig
end

##
SqsGFMC = fetch.([Threads.@spawn getSq(res,150÷nBra,nThermal) for res in results])
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
#___________Periodic Boundaries_______________________
function getPeriodic(parent)
    state = parent |> Array
    SW.SpinConfig(SW.PeriodicMatrix(state), parent.S)
end

S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
S_ED = getPeriodic(SW.getStairCase(size(S,1)))
# S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2)
# S_ED = SW.getStairCase(size(S,1))
HStair = SW.generateHilbertSpace(S_ED)

ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1]
HConfs = SW.spinConfig.(HStair.AllStates,Ref(S_ED),Ref(HStair.plaqMapping))
magEx = SW.getMagnetization(HConfs, v0)
SqEx = SW.getStructureFac(HConfs,v0)
##
nThermal = 100
nBra = 2
ψG = SW.PlaquetteNumberGuidingFunction(0.197)
# ψG(N) = 1
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,4,(nThermal+1300)÷nBra,nBra,ψG,1) for _ in 1:6])

##
plotEnergies(results,nBra,nThermal,E0,Emin=E0-1e-1,Emax=E0+2e-1)

##___________ StraightForwardWalking _______________________
allPlaqs = SW.getApplicablePlaquettes(S_ED)
refPlaq = allPlaqs[begin + 13]
# pairPlaqs = filter(!=(refPlaq),allPlaqs)
pairPlaqs = copy(allPlaqs)
exactCorrs = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,refPlaq,Pj) for Pj in pairPlaqs]
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    points = Point.(pairPlaqs)
    sizefunc(x) = x*100
    scatter!(ax,Point(refPlaq), marker = '×',markersize = 60, color = :red)
    scatter!(ax,points, markersize = sizefunc.(exactCorrs),colormap = :viridis, color = exactCorrs)
    fig
end
##
resSFW = fetch.([Threads.@spawn SW.measure_B2_correlators(S,res.SaveConfigs,5,nBra,refPlaq,ψG,1) for res in results])
##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,20) for res in results]

GFMCPlaqs = collect(SW.plaquetteIterator(S))

obs = [[SW.get_observables_sfw(Gnp,res[:,:,j],mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resSFW,results) ]
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(GFMCPlaqs)

    corrEnd = last.(mean(obs))

    localCorr = only(findfirst(==(refPlaq),GFMCPlaqs))
    corrEnd[localCorr] /= 2

    exactCorrsRescale = copy(exactCorrs)
    localCorr = only(findfirst(==(refPlaq),pairPlaqs))
    exactCorrsRescale[localCorr] /= 2

    sizefunc(x) = x*30*7
    # scatter!(ax,Point(refPlaq), marker = '×',markersize = 60, color = :red)
    
    points = Point.(pairPlaqs)
    sizefunc2(x) = sizefunc(x)
    scatter!(ax,points, markersize = sizefunc2.(exactCorrsRescale),colormap = :viridis, color = sizefunc2.(exactCorrsRescale),marker = '●',alpha = 1)
    
    scatter!(ax,pointsGFMC, markersize = sizefunc.(corrEnd),colormap = :viridis, color = :red,alpha = 0.0,marker = '●',strokewidth = 0.8000,strokecolor = :red) 
    # scatter!(ax,points, markersize = sizefunc2.(exactCorrs),colormap = :viridis, color = exactCorrs,marker = '∘')
    # scatter!(ax,points, markersize = sizefunc2.(exactCorrs),colormap = :viridis, color = exactCorrs,alpha = 0.4)
    fig
end
##
Plaq2 = (8,7)
# Plaq2 = refPlaq

gfmcPlaq = only(findfirst(==(Plaq2),GFMCPlaqs))
obsArr = stack(stack(obs))[:,gfmcPlaq,:]
errorbars(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],sqrt.(var(obsArr,dims=2))[:])
lines!(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],linewidth = 0.5)
exactCorr = exactCorrs[only(findfirst(==(Plaq2),pairPlaqs))]
hlines!([exactCorr],color = :red)
current_figure()
