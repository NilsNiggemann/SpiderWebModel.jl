import Pkg
Pkg.activate("Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

function plotEnergies(results,nBra,nThermal)
    ens = [SW.getEnergies(res.TotalWeights[nThermal:end],res.energies[nThermal:end],1,250÷nBra) for res in results]

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
        # hlines!([E0],color = :red,label = L"exact $$")
        axislegend(ax,merge=true)
        # xlims!(ax,0.5,last(proj))
        # ylims!(ax,Emin,Emax)
        # save("Application/exactFig/GFMCEnergy.png",fig)
        fig
    end
end

##
#___________Periodic Boundaries_______________________
function getPeriodic(parent)
    state = parent |> Array
    SW.SpinConfig(SW.PeriodicMatrix(state), parent.S)
end

S = SW.stencilConfig(parent(SW.getStairCase(12)),1/2;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
S_ED = getPeriodic(SW.getStairCase(size(S,1)))
# S_ED = SW.getStairCase(size(S,1))
# S = SW.stencilConfig(parent(SW.getStairCase(7)),1/2)
# S_ED = SW.getStairCase(size(S,1))
HStair = SW.generateHilbertSpace(S_ED)

ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1]
HConfs = getPeriodic.(SW.spinConfig.(HStair.AllStates,Ref(S_ED),Ref(HStair.plaqMapping)))
magEx = SW.getMagnetization(HConfs, v0)
SqEx = SW.getStructureFac(HConfs,v0)
##
nThermal = 100
nBra = 2
ψG = SW.PlaquetteNumberGuidingFunction(0.197)
# ψG(N) = 1
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,4,(nThermal+4800)÷nBra,nBra,ψG,1) for _ in 1:6])

##
plotEnergies(results,nBra,nThermal)


##___________ StraightForwardWalking _______________________
# allPlaqs = collect(SW.plaquetteIterator(S))
allPlaqs = SW.getApplicablePlaquettes(S_ED)
refPlaq = first(SW.getApplicablePlaquettes(S_ED))
# pairPlaqs = filter(!=(refPlaq),allPlaqs)
pairPlaqs = copy(allPlaqs)
exactBB = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,refPlaq,Pj) for Pj in pairPlaqs]
exactB = [SW.getBi_square(HStair.AllStates,HStair.plaqMapping,v0,Pj) for Pj in pairPlaqs]
refPlaqPos = only(findfirst(==(refPlaq),pairPlaqs))
exactCorrs = exactBB .- exactB[refPlaqPos] .*exactB 
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
BOp = SW.PlaquetteFlipOperator(S)
resB = fetch.([Threads.@spawn SW.measure_operator(S,res.SaveConfigs,8,nBra,BOp,ψG,1) for res in results])
##
BBOp = SW.BBOperator(S,refPlaq)
resBB = fetch.([Threads.@spawn SW.measure_operator(S,res.SaveConfigs,8,nBra,BBOp,ψG,1) for res in results])

##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,10) for res in results]

GFMCPlaqs = collect(SW.plaquetteIterator(S))

BVals = [[SW.get_observables_sfw(Gnp,res[:,:,j],mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resB,results) ]
BBVals = [[SW.get_observables_sfw(Gnp,res[:,:,j],mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resBB,results) ]
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(GFMCPlaqs)

    BBEnd = last.(mean(BBVals)) 
    BEnd = last.(mean(BVals))
    localCorr = only(findfirst(==(refPlaq),GFMCPlaqs))
    corrEnd = BBEnd .- BEnd[localCorr] .*BEnd
    
    corrEnd[localCorr] /= 2

    exactCorrsRescale = copy(exactCorrs)
    localCorr = only(findfirst(==(refPlaq),pairPlaqs))
    exactCorrsRescale[localCorr] /= 2

    sizefunc(x) = x*30*25
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
Plaq2 = (3,4)
# Plaq2 = rand(pairPlaqs)
# Plaq2 = refPlaq

gfmcPlaq = only(findfirst(==(Plaq2),GFMCPlaqs))
obsArr = stack(stack(BBVals))[:,gfmcPlaq,:]
errorbars(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],sqrt.(var(obsArr,dims=2))[:])
lines!(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],linewidth = 0.5)
exactCorr = exactCorrs[only(findfirst(==(Plaq2),pairPlaqs))]
hlines!([exactCorr],color = :red)
# ylims!(exactCorr -5e-1,exactCorr + 5e-1)
current_figure()
##
#___________Spin-1_______________________
ψG = SW.PlaquetteNumberGuidingFunction(0.15)
S = SW.stencilConfig(zeros(8,8),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
##
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,8,52000÷nBra,nBra,ψG,1) for _ in 1:6])
##
plotEnergies(results,nBra,nThermal)
##
refPlaq = (5,6)
wrap_idx(x,L) = abs(x) >= L ÷ 2 ? x - sign(x)*L : x
# function wrap_idx(x,L)
#     if x <= -(L ÷2)
#         return x + L
#     elseif x > (L ÷2)
#         return x - L
#     end
#     return x
# end
function isReducedPlaq(P,refPlaq,L)
    x,y = P .- refPlaq
    x = wrap_idx(x,L)
    y = wrap_idx(y,L)
    return x>=0 && y>=0 && x>=y
end
allPlaqs = collect(SW.plaquetteIterator(S))
reducedPlaqs = filter(P->isReducedPlaq(P,refPlaq,size(S,1)),allPlaqs)
# scatter(Point.(allPlaqs))
# scatter!(Point.(reducedPlaqs),color = :red)
# current_figure()
##
BOp = SW.PlaquetteFlipOperator(S)
resB = fetch.([Threads.@spawn SW.measure_operator(S,res.SaveConfigs,50,nBra,BOp,ψG,1) for res in results])
##
BBOp = SW.BBOperator(S,refPlaq)
resBB = fetch.([Threads.@spawn SW.measure_operator(S,res.SaveConfigs,50,nBra,BBOp,ψG,1) for res in results])

##

Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,50) for res in results]

BBVals = [[SW.get_observables_sfw(Gnp,res[:,:,j],mean(result.TotalWeights)) for j in eachindex(allPlaqs)] for (Gnp,res,result) in zip(Gnps,resBB,results) ]
BVals = [[SW.get_observables_sfw(Gnp,res[:,:,j],mean(result.TotalWeights)) for j in eachindex(allPlaqs)] for (Gnp,res,result) in zip(Gnps,resB,results) ]

##
using SpiderWebModel.StaticArrays


function mapPlaqToReduced(P,refPlaq,L)
    x´ = abs(P[1] - refPlaq[1])
    y´ = abs(P[2] - refPlaq[2])
    wrap(x) = wrap_idx(x,L)
    if x´ > y´
        return wrap.((x´,y´)) .+refPlaq
    else
        return wrap.((y´,x´)) .+refPlaq
    end
end

function getBBCorrelator(BBVals,BVals,redPlaqs,allPlaqs,refPlaq)
    nth(x) = x[30]
    BBEnd = nth.(mean(BBVals)) 
    BEnd =  nth.(mean(BVals))
    localCorr = only(findfirst(==(refPlaq),reducedPlaqs))
    BBCorrelatorRaw = BBEnd .- BEnd[localCorr] .*BEnd
    return BBCorrelatorRaw = BBEnd .- BEnd[localCorr] .*BEnd
    # L = maximum(maximum,allPlaqs)
    # BBCorrelator = empty(BBCorrelatorRaw)
    # for P in allPlaqs
    #     P´ = mapPlaqToReduced(P,refPlaq,L)
    #     i = findfirst(==(P´),redPlaqs)
    #     if i === nothing
    #         # error("Plaquette $P not found")
    #         push!(BBCorrelator,0)
    #     else
    #         push!(BBCorrelator,BBCorrelatorRaw[i])
    #     end
    # end
    # return BBCorrelator
end
function getBBCorrelator(BBVals,BVals,refPlaq)
    nth(x) = x[30]
    BBEnd = nth.(mean(BBVals)) 
    BEnd =  nth.(mean(BVals))
    localCorr = only(findfirst(==(refPlaq),reducedPlaqs))
    BBCorrelatorRaw = BBEnd .- BEnd[localCorr] .*BEnd
    return BBCorrelatorRaw = BBEnd .- BEnd[localCorr] .*BEnd
end

# inds = findall(P->isReducedPlaq(P,refPlaq,size(S,1)),allPlaqs)
# BBVals = [x[inds] for x in BVals]
# BVals = [x[inds] for x in BVals]

BBCorrelator = getBBCorrelator(BBVals,BVals,refPlaq)
    
# BBCorrelator = getBBCorrelator(BBVals,BVals,reducedPlaqs,allPlaqs,refPlaq)
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(allPlaqs)

    corrEnd = copy(BBCorrelator)
    # corrEnd = last.(mean(BBVals))
    perm = sortperm(pointsGFMC .- Point(refPlaq),rev =false, by = SW.norm)
    corrEnd = corrEnd[perm]
    pointsGFMC = pointsGFMC[perm]
    # corrEnd = copy( BEnd[localCorr] .*BEnd)
    # corrEnd = copy( BBEnd)
    sizefunc(x) = abs(x)*30*8
    sizes = sizefunc.(corrEnd) 

    corrEnd[localCorr] /= 1
    
    points = Point.(pairPlaqs)
    markerfunc(x) = x>0 ? '●' : '○'
    scatter!(ax,pointsGFMC, markersize = sizes,colormap = :viridis, color = sizefunc.(corrEnd),alpha = 1.0,marker = markerfunc.(corrEnd))
    fig
end
##
with_theme(theme_PiTicks()) do
    
    obsMat = zeros(size(S))
    obsVals = BBCorrelator
    for (i,I) in enumerate(SW.plaquetteIterator(S))
        obsMat[I...] = obsVals[i]
    end
    FT = SW.LatticeFFTs.fft(obsMat)
    FT[1,1] = NaN
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    k = 0:size(FT,1) ./ size(FT,1) .*2pi
    heatmap!(ax,k,k,real(FT))
    fig
end
##


function FTPlaq(rPlaq,Vals,k)
    res = 0.
    for (r,Val) in zip(rPlaq,Vals)
        res += Val * cos(k'*r)
    end
    res
end

with_theme(theme_PiTicks()) do
    Tinv = SW.SA[
        1 -1;
        1 1
    ]
    ri = [SW.SVector(r .- refPlaq) for r in SW.plaquetteIterator(S)]
    rPlaq = [mapToPlaquetteBasis(r) for r in ri]

    # return scatter(Point.(ri),color = BBCorrelator,markersize = 150 * abs.(BBCorrelator))
    k = LinRange(-pi,pi,500)

    FT = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    hm = heatmap!(ax,k,k,FT)
    Colorbar(fig[1,2],hm)
    fig
end
##
function mapToPlaquetteBasis(I)
    r = SW.SA[I...]
    T = SW.SA[
        1 -1;
        1  1
    ]
    return T*r .÷2
end
with_theme(theme_PiTicks()) do
    Tinv = SW.SA[
        1 -1;
        1 1
    ]
    # ri = [SW.SVector(r .- refPlaq) for r in SW.plaquetteIterator(S)]
    ri = [SW.SVector(r) for r in SW.plaquetteIterator(S)]
    rPlaq = [mapToPlaquetteBasis(r) .+ Point(refPlaq) for r in ri]
    # return sort(rPlaq)
    obsMat = zeros(size(S) .+10)
    for (i,I) in enumerate(rPlaq)
        ii,jj = I
        obsMat[ii+5,jj+5 ] = BBCorrelator[i]
    end
    return heatmap(obsMat)

    FT = SW.LatticeFFTs.fft(obsMat)
    FT[1,1] = NaN
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    k = 0:size(FT,1) ./ size(FT,1) .*2pi
    heatmap!(ax,k,k,real(FT))
    fig
end
