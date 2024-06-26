import Pkg
# Pkg.activate("../Application/")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
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
GC.gc()
GC.gc()

outfile = "Data/temp/S12_8/"
rm(outfile;recursive=true,force=true)
mkpath(outfile)
nThermal = 1000
nBra = 3
ψG = SW.fullVariationalFunction(S,0.197)

# ψG(N) = 1
# DT = SW.ContinuousTimeMethod(0.1,3,-E0)
DT = SW.DiscreteTimeMethod(0.,3,-E0)
# stochReconfRes = SW.stochastic_reconfiguration(S,DT,i->round(Int,1000+ 200*i),ψG,50,0.6,SW.IterativeSRSolver();Nwalkers = 6*8,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000)
# ψG = typeof(ψG)(stochReconfRes.params)
# DT = SW.DiscreteTimeMethod(0,3,E0)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,20,2000,ψG,equilibration_steps=100,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:6])

##
plotEnergies(results,DT,E0,nThermal=10,normalize=true,Emin=NaN)

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
resB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,5,BOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])
# resB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,5,BOp,ψG,1;outfile = string(outfile,i,".h5")) for (i,res) in enumerate(results)])
##
BBOp = SW.BBOperator(S,refPlaq)
resBB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,5,BBOp,ψG,SW.getApplicablePlaquettes(S)) for (i,res) in enumerate(results)])

##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,5) for res in results]

GFMCPlaqs = SW.getApplicablePlaquettes(S)
# GFMCPlaqs = collect(SW.plaquetteIterator(S))

BVals = [[SW.get_observables_sfw(Gnp,res[:,j,:]',mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)[1:1]] for (Gnp,res,result) in zip(Gnps,resB,results) ]
BBVals = [[SW.get_observables_sfw(Gnp,res[:,j,:]',mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resBB,results) ]
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(GFMCPlaqs)
    nth(x) = x[end]
    BBEnd = nth.(mean(BBVals)) 
    BEnd = only(nth.(mean(BVals)))
    localCorr = only(findfirst(==(refPlaq),GFMCPlaqs))
    corrEnd = BBEnd .- BEnd^2
    
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
Plaq2 = (4,7)
# Plaq2 = rand(pairPlaqs)
# Plaq2 = refPlaq

gfmcPlaq = only(findfirst(==(Plaq2),GFMCPlaqs))
obsArr = stack(stack(BBVals))[:,gfmcPlaq,:] 
errorbars(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],sqrt.(var(obsArr,dims=2))[:])
lines!(eachindex(nBra .* obsArr[:,1]),mean(obsArr,dims=2)[:],linewidth = 0.5)
exactCorr = exactBB[only(findfirst(==(Plaq2),pairPlaqs))]
hlines!([exactCorr],color = :red)
# ylims!(exactCorr -5e-1,exactCorr + 5e-1)
current_figure()
##
#___________Spin-1_______________________
S = SW.stencilConfig(zeros(16,16),1;boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional())
DT = SW.DiscreteTimeMethod(0.,2,prod(size(S)))
ψG = SW.fullVariationalFunction(S,0.15)
stochReconfRes = SW.stochastic_reconfiguration(S,DT,i->round(Int,200+ 5*i),ψG,30,0.8,SW.IterativeSRSolver();Nwalkers = 6*20,rel_tolerance=1e-8,equilibration_steps=100,pre_equilibration_steps=10_000)
##
ψG = typeof(ψG)(stochReconfRes.params)
DT = SW.DiscreteTimeMethod(0.,8,0.2658*prod(size(S)))

@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,DT,20,1000,ψG,equilibration_steps=5000,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:6])
##
plotEnergies(results,nBra,p=100;normalize=true)
##
# scatter(Point.(allPlaqs))
# scatter!(Point.(reducedPlaqs),color = :red)
# current_figure()
refPlaq = SW.getCentralPlaquette(S)
symReduc = SW.symmetryReducePlaquettes(S,refPlaq)
GFMCPlaqs = collect(SW.plaquetteIterator(S))[symReduc.uniqueInds]
allPlaqs = collect(SW.plaquetteIterator(S))

##
outfile = "Data/temp/S1/L=$(size(S,1))"
mkpath(outfile)
BOp = SW.PlaquetteFlipOperator(S)
resB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,1,BOp,ψG,[refPlaq];outfile = outfile*"$i.h5") for (i,res) in enumerate(results)])
##
BBOp = SW.BBOperator(S,refPlaq)
# resBB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,1,BBOp,ψG,GFMCPlaqs) for (i,res) in enumerate(results)])
resBB = fetch.([Threads.@spawn SW.measure_operator(S,DT,res.SaveConfigs,1,BBOp,ψG,GFMCPlaqs;outfile = outfile*"$i.h5") for (i,res) in enumerate(results)])


##
Gnps = [SW.precomputeNormalizedAccWeight(res.TotalWeights,1,1) for res in results]

BBVals = [[SW.get_observables_sfw(Gnp,res[:,j,:]',mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resBB,results) ]
BVals = [SW.get_observables_sfw(Gnp,res[:,begin,:]',mean(result.TotalWeights)) for (Gnp,res,result) in zip(Gnps,resB,results) ]

##
using SpiderWebModel.StaticArrays

function getBBCorrelator(BBVals,BVals,symReduc;index = lastindex(BBVals[1][1]))

    BBCorrelatorRaw = getBBCorrelator(BBVals,BVals,index)

    BBCorrelator = similar(BBCorrelatorRaw, length(symReduc.indicesMapping))
    for (i,k) in enumerate(symReduc.indicesMapping)
        BBCorrelator[i] = BBCorrelatorRaw[k]
    end
    return BBCorrelator
end

function getBBCorrelator(BBVals,BVals,index::Int)
    nth(x) = x[index]
    BBEnd = nth.(mean(BBVals))
    BEnd =  nth(mean(BVals))
    BBCorrelatorRaw = BBEnd .- BEnd^2
    return BBCorrelatorRaw
end

# inds = findall(P->isReducedPlaq(P,refPlaq,size(S,1)),allPlaqs)
# BBVals = [x[inds] for x in BVals]
# BVals = [x[inds] for x in BVals]
BBCorrelator = getBBCorrelator(BBVals,BVals,symReduc,index = 1)
# BBCorrelator = getBBCorrelator(BBVals,BVals,1)

# BBCorrelator = getBBCorrelator(BBVals,BVals,reducedPlaqs,allPlaqs,refPlaq)
##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(collect(SW.plaquetteIterator(S)))

    corrEnd = copy(BBCorrelator)
    # corrEnd = last.(mean(BBVals))
    perm = sortperm(pointsGFMC .- Point(refPlaq),rev =false, by = SW.norm)
    corrEnd = corrEnd[perm]
    pointsGFMC = pointsGFMC[perm]
    # corrEnd = copy( BEnd[localCorr] .*BEnd)
    # corrEnd = copy( BBEnd)
    sizefunc(x) = abs(x)*30*10
    localScale = 2
    corrEnd[1] /= localScale
    sizes = sizefunc.(corrEnd) 
    corrEnd[1] *= localScale
    
    points = Point.(pairPlaqs)
    markerfunc(x) = x>0 ? '●' : '○'
    scatter!(ax,pointsGFMC, markersize = sizes,colormap = :viridis, color = sizefunc.(corrEnd),alpha = 1.0,marker = markerfunc.(corrEnd))
    fig
end


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
    ri = [T * SW.SVector(r .- refPlaq) for r in SW.plaquetteIterator(S)]
    # rPlaq = [mapToPlaquetteBasis(r) for r in ri]
    capfilter(x) = min(abs(x),30)

    # return scatter(Point.(ri),color = capfilter.(200 .*BBCorrelator./maximum(BBCorrelator)),markersize = capfilter.(400 * BBCorrelator./maximum(BBCorrelator)),axis =(;aspect=1,xticks=SimpleTicks(-6:6), yticks=SimpleTicks(-6:6)))

    # return scatter(Point.(ri),color = BBCorrelator,markersize = 150 * abs.(BBCorrelator))
    k = LinRange(-pi,pi,200)

    FT = [FTPlaq(ri,BBCorrelator,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
    fig = Figure(size = 1.3 .*(500,250))
    ax = Axis(fig[1,1];aspect = 1)
    ax2 = Axis(fig[1,2];aspect = 1)

    hm = heatmap!(ax,k,k,FT)

    photw = ω_photon.(Iterators.product(k,k))
    photw .*= maximum(FT)/maximum(photw)

    heatmap!(ax2,k,k,photw,colormap = :viridis,colorrange = extrema(FT))
    Colorbar(fig[1,3],hm)
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

    FT = SW.FFTW.fft(obsMat)
    FT[1,1] = NaN
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    k = 0:size(FT,1) ./ size(FT,1) .*2pi
    heatmap!(ax,k,k,real(FT))
    fig
end
