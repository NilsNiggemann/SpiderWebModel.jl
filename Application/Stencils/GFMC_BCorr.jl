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
function getPeriodic(mat,Spin)
    SW.SpinConfig(SW.PeriodicMatrix(0.5*parent(mat)), Spin)
end

mu = 0.3


S = SW.stencilConfig(parent(SW.periodicState6x6_3(12)),1/2;boundaryCondition = :periodic)
# SW.flipSpinsAlongDiagonal!(S,6,1)
# SW.flipSpinsAlongRow!(S,3,1)

S_ED = getPeriodic(S,1/2)
SW.plotApplPlaquettes(S)
##
# S_ED = SW.getStairCase(size(S,1))
# S = SW.stencilConfig(parent(SW.getStairCase(7)),1/2)
# S_ED = SW.getStairCase(size(S,1))
HStair = SW.generateHilbertSpace(S_ED)
SW.addRKPotential!(HStair,mu)
ExSol = SW.SolveHKrylov(HStair.H)
E0 = ExSol.values[1]
v0 = ExSol.vectors[1]
HConfs = getPeriodic.(SW.spinConfig.(HStair.AllStates,Ref(S_ED),Ref(HStair.plaqMapping)),1/2)
magEx = SW.getMagnetization(HConfs, v0)
SqEx = SW.getStructureFac(HConfs,v0)
##

CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = E0,Hxx = SW.Hxx_RK(mu))
CTSR = SW.ContinuousTimeMethod(30*CT.τ,w_avg_estimate = CT.w_avg_estimate,Hxx = CT.Hxx)

##

nThermal = 1000
nBra = 3
ψG = SW.SimpleJastrowFunction(S)
SW.rand!(SW.getNonSymmetric(ψG),1e-3)
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,20,ψG,100,3e-3,SW.IterativeSRSolver();Nwalkers = 20,reconfigure = false,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=10)
SW.get_params(ψG) .= stochReconfRes.params
plotVarEn(stochReconfRes)
##
# ψG(N) = 1
# CT = SW.ContinuousTimeMethod(0.1,3,-E0)
# stochReconfRes = SW.stochastic_reconfiguration(S,CT,i->round(Int,1000+ 200*i),ψG,50,0.6,SW.IterativeSRSolver();Nwalkers = 6*8,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000)
# ψG = typeof(ψG)(stochReconfRes.params)
# CT = SW.DiscreteTimeMethod(0,3,E0)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,20*1,300,ψG,equilibration_steps=100,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:20])

##
plotEnergies(results,CT,E0,nThermal=10,normalize=true,Emin = E0-1e-2,Emax = E0+1.5e-2)

##___________ StraightForwardWalking _______________________
# allPlaqs = collect(SW.plaquetteIterator(S))
allPlaqs = SW.getApplicablePlaquettes(S_ED)
refPlaq = first(allPlaqs)
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
# resB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,10,BOp,ψG) for (i,res) in enumerate(results)])
# BVals = stack(stack([[SW.get_observables_sfw(Gnp,res[:,j,:]',mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)[1:1]] for (Gnp,res,result) in zip(Gnps,resB,results) ]))
BVals = stack(fetch.([Threads.@spawn SW.measureObservables(S,BOp,allPlaqs,5,20*1,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:10]))
##
BOp_rand = SW.RandomPlaquetteFlipOperator(S)
BValsRand = stack([SW.measureObservables(S,BOp_rand,[nothing],5,20*3,3000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:10]) ./ length(allPlaqs)
##
# BVals = fetch.([Threads.@spawn SW.measureObservables(S,BOp,allPlaqs,5,20*1,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:20])
##
BBOp = SW.BBOperator(S,refPlaq)
# GFMCPlaqs = SW.getApplicablePlaquettes(S)
GFMCPlaqs = allPlaqs

# resBB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,3,BBOp,ψG,SW.getApplicablePlaquettes(S)) for (i,res) in enumerate(results)])
# BBVals = stack(stack([[SW.get_observables_sfw(Gnp,res[:,j,:]',mean(result.TotalWeights)) for j in eachindex(GFMCPlaqs)] for (Gnp,res,result) in zip(Gnps,resBB,results) ]))

BBVals = stack([SW.measureObservables(S,BBOp,GFMCPlaqs,5,20*2,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:10]) 
##

# GFMCPlaqs = collect(SW.plaquetteIterator(S))

##
let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(GFMCPlaqs)
    nth(x) = x[end]
    BBEnd = dropmean(BBVals,dims=3)[end,:]
    BEnd = dropmean(BVals,dims=3)[end,:]
    localCorr = only(findfirst(==(refPlaq),GFMCPlaqs))
    # corrEnd = BBEnd .- BEnd[localCorr].*BEnd
    corrEnd = BEnd
    # corrEnd = BEnd^2
    
    corrEnd[localCorr] /= 2

    # exactCorrsRescale = copy(exactCorrs)
    exactCorrsRescale = copy(exactB)

    # exactCorrsRescale = exactB
    localCorr = only(findfirst(==(refPlaq),pairPlaqs))
    exactCorrsRescale[localCorr] /= 2
    
    sizefunc(x::AbstractArray) = x.*100 ./ maximum(x)
    sizefunc2(x) = sizefunc(x)
    # scatter!(ax,Point(refPlaq), marker = '×',markersize = 60, color = :red)
    
    points = Point.(pairPlaqs)
    scatter!(ax,points, markersize = sizefunc(exactCorrsRescale),colormap = :viridis, color = sizefunc(exactCorrsRescale),marker = '●',alpha = 1)
    
    scatter!(ax,pointsGFMC, markersize = sizefunc(corrEnd),colormap = :viridis, color = :red,alpha = 0.0,marker = '●',strokewidth = 0.8000,strokecolor = :red) 
    # scatter!(ax,points, markersize = sizefunc2.(exactCorrs),colormap = :viridis, color = exactCorrs,marker = '∘')
    # scatter!(ax,points, markersize = sizefunc2.(exactCorrs),colormap = :viridis, color = exactCorrs,alpha = 0.4)
    fig
end
##


Plaq2 = (4,7)
# Plaq2 = refPlaq
with_theme(theme_SimpleTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1])
    for Plaq2 in pairPlaqs
        gfmcPlaq = only(findfirst(==(Plaq2),GFMCPlaqs))
        # obsArr = BBVals[:,gfmcPlaq,:] .- BVals[:,1,:].^2
        obsArr = BVals[:,1,:].^2

        # exactCorr = exactBB[only(findfirst(==(Plaq2),pairPlaqs))] - exactB[refPlaqPos]^2
        exactCorr = exactB[refPlaqPos]^2

        errorbars!(eachindex(nBra .* obsArr[:,1]),dropmean(obsArr,dims=2),dropstd(obsArr,dims=2))
        l = lines!(eachindex(nBra .* obsArr[:,1]),dropmean(obsArr,dims=2),linewidth = 0.5)

        # Ex_diff = dropmean(obsArr,dims=2) .- exactCorr
        # Ex_diff_std = dropstd(obsArr,dims=2)
        # l = errlines!(ax,eachindex(Ex_diff),Ex_diff,Ex_diff_std,linewidth = 2)
        hlines!([exactCorr],color = (l.color[],0.5),linewidth = 1,linestyle = :dash)

    end
    # ylims!(exactCorr -5e-1,exactCorr + 5e-1)
    current_figure()
end
##

function PlaqSumFT(S,BBCorr,refPlaq,allPlaqs,shift = 0)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                rx,ry = P .- refPlaq
                BBq += BBCorr[iP] * cos(kx*rx + ky*ry+shift)
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

function PlaqSumFTPairs(S,BBCorr::AbstractMatrix,allPlaqs,shift = 0)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                for (jP,P2) in enumerate(allPlaqs)
                    rx,ry = P .- P2
                    BBq += BBCorr[iP,jP] * cos(kx*rx + ky*ry+shift)
                end
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

##
BBQOp = SW.BBqOperator()
qvals = let 
    qx = trueMomenta(0,2pi,size(S,1))
    qy = trueMomenta(0,2pi,size(S,2))
    qvals = [SA[qx,qy] for (qx,qy) in Iterators.product(qx,qy)]
end
BBQOp = SW.BBqOperator()
##
BBQ = stack(fetch.([Threads.@spawn SW.measureObservables(S,BBQOp,qvals,5,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))

##
exactBB_pairs = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,Pi,Pj) for Pi in allPlaqs, Pj in allPlaqs] ./ length(allPlaqs)

with_theme(theme_PiTicks()) do 
    fig = Figure(size = (400,600))
    ax1 = Axis(fig[1,1];aspect = 1,title = "exact BBq")

    ax2 = Axis(fig[2,1];aspect = 1,title = "GFMC BBq op")
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    FTex = PlaqSumFTPairs(S,exactBB_pairs,allPlaqs)
    hm = heatmap!(ax1,kx,ky,FTex,colormap = :viridis)
    Colorbar(fig[1,2],hm)
    fig

    FTmean = dropmean(BBQ,dims=3)[end,:] ./ length(allPlaqs)
    FT = zeros(length(kx),length(ky))
    FT[:] .= FTmean
    FT .-= FT[1,1]*0.5
    hm = heatmap!(ax2,kx,ky,FT,colormap = :viridis)
    Colorbar(fig[2,2],hm)
    fig
end
##___________Bq Operator test_______________________

BQOp = SW.BqOperator()
BQ = stack(fetch.([Threads.@spawn SW.measureObservables(S,BQOp,qvals,1,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))

##
BQSin = stack(fetch.([Threads.@spawn SW.measureObservables(S,BQOp,qvals,1,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))
##
with_theme(theme_PiTicks()) do 
    fig = Figure(size = (400,600))
    ax1 = Axis(fig[1,1];aspect = 1,title = "exact Bq")

    ax2 = Axis(fig[2,1];aspect = 1,title = "GFMC Bq op")
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    FTex = PlaqSumFT(S,exactB,SW.getCentralPlaquette(S),allPlaqs)
    hm = heatmap!(ax1,kx,ky,FTex,colormap = :viridis)
    Colorbar(fig[1,2],hm)
    fig

    FTmean = dropmean(BQ,dims=3)[end,:] #./ length(allPlaqs)
    FT = zeros(length(kx),length(ky))
    FT[:] .= FTmean
    FT .-= FT[1,1]*0.5
    hm = heatmap!(ax2,kx,ky,FT,colormap = :viridis)
    Colorbar(fig[2,2],hm)
    fig
end
## _________together_______________________
BiBj = [exactB[i]*exactB[j] for i in eachindex(exactB), j in eachindex(exactB)] ./length(allPlaqs)
# exactBB_pairs_decor = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,Pi,Pj) for Pi in allPlaqs, Pj in allPlaqs] ./ length(allPlaqs)


with_theme(theme_PiTicks()) do 
    fig = Figure(size = (400,600))
    ax1 = Axis(fig[1,1];aspect = 1,title = "exact B(q)B(-q)")

    ax2 = Axis(fig[2,1];aspect = 1,title = "GFMC Bq op")
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    FTex = PlaqSumFTPairs(S,exactBB_pairs .- BiBj ,allPlaqs)
    # FTex = PlaqSumFTPairs(S,exactBB_pairs ,allPlaqs) .-( PlaqSumFT(S,exactB,SW.getCentralPlaquette(S),allPlaqs)./length(allPlaqs)).^2


    hm = heatmap!(ax1,kx,ky,FTex,colormap = :viridis)
    Colorbar(fig[1,2],hm)

    FTmean = dropmean(BQ,dims=3)[end,:] #./ length(allPlaqs)
    FT = zeros(length(kx),length(ky))
    FT[:] .= FTmean
    FT .-= FT[1,1]*0.5
    hm = heatmap!(ax2,kx,ky,FT,colormap = :viridis)
    Colorbar(fig[2,2],hm)
    fig
end

##
#___________Spin-1_______________________
S = SW.stencilConfig(zeros(8,8),1;boundaryCondition = :periodic)
# S .= SW.periodicStateDenseLoops(size(S,1))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.5))
CTSR = SW.ContinuousTimeMethod(10*CT.τ,Hxx = CT.Hxx)

ψG = SW.SimpleJastrowFunction(S)
# ψG = SW.RKFunction()
##
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,10,ψG,300,1e-3,SW.IterativeSRSolver();Nwalkers = 2*20,reconfigure = false,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=10_000)
SW.get_params(ψG) .= stochReconfRes.params
plotVarEn(stochReconfRes,movavg = 30)
##
# scatter(Point.(allPlaqs))
# scatter!(Point.(reducedPlaqs),color = :red)
# current_figure()
refPlaq = SW.getCentralPlaquette(S)
symReduc = SW.symmetryReducePlaquettes(S,refPlaq)
# GFMCPlaqs = collect(SW.plaquetteIterator(S))[symReduc.uniqueInds]
GFMCPlaqs = collect(SW.plaquetteIterator(S))
allPlaqs = collect(SW.plaquetteIterator(S))

##
BOp = SW.PlaquetteFlipOperator(S)
# resB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,10,BOp,ψG,[refPlaq]) for (i,res) in enumerate(results)])
BVals = fetch.([Threads.@spawn SW.measureObservables(S,BOp,GFMCPlaqs,10,20*3,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:15])

##
BBOp = SW.BBOperator(S,refPlaq)
BBVals = fetch.([Threads.@spawn SW.measureObservables(S,BBOp,GFMCPlaqs,10,20*3,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:15])
##
BBOp2 = SW.BBOperator(S,refPlaq.+ (-1,1))
BBVals2 = fetch.([Threads.@spawn SW.measureObservables(S,BBOp,GFMCPlaqs,10,20*3,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:15])

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

function getBBCorrelator(BBVals,BVals,refPlaq,allPlaqs,index::Int)
    BBEnd = BBVals[index,:]
    refPlaqInd = only(findfirst(==(refPlaq),allPlaqs))
    BEnd =  BVals[index,:]
    BBCorrelatorRaw = BBEnd .- BEnd.*BEnd[refPlaqInd]
    return BBCorrelatorRaw
end

# inds = findall(P->isReducedPlaq(P,refPlaq,size(S,1)),allPlaqs)
# BBVals = [x[inds] for x in BVals]
# BVals = [x[inds] for x in BVals]
# BBCorrelator = getBBCorrelator(BBVals,BVals,symReduc,index = 1)

BBCorrelator = getBBCorrelator.(0BBVals, BVals,Ref(refPlaq),Ref(allPlaqs),3)
# BBCorrelator2 = getBBCorrelator.(0 .*BBVals2, BVals,Ref(refPlaq .+ (-1,1)),Ref(allPlaqs),3)

# BBCorrelator = getBBCorrelator(BBVals,BVals,reducedPlaqs,allPlaqs,refPlaq)
##
function FTPlaq(rPlaq,Vals,k)
    res = 0.
    for (r,Val) in zip(rPlaq,Vals)
        res += Val * cos(k'*r)
    end
    res
end

let 
    fig = Figure()
    ax = Axis(fig[1,1];SW.getConfigAxis(S)...,backgroundcolor = :white)
    pointsGFMC = Point.(collect(SW.plaquetteIterator(S)))

    corrEnd = mean(BBCorrelator)
    # corrEnd = dropmean(BBCorrRK,dims=2)
    # corrEnd = last.(mean(BBVals))
    perm = sortperm(pointsGFMC .- Point(refPlaq),rev =false, by = SW.norm)
    corrEnd = corrEnd[perm]
    pointsGFMC = pointsGFMC[perm]
    # corrEnd = copy( BEnd[localCorr] .*BEnd)
    # corrEnd = copy( BBEnd)

    localScale = 5

    corrEnd[1] /= localScale
    
    TotNorm = maximum(abs,corrEnd)

    sizefunc(x) = abs(x)*30*2 / TotNorm
    sizes = sizefunc.(corrEnd) 
    
    corrEnd[1] *= localScale

    markerfunc(x) = x>0 ? '●' : '○'
    scatter!(ax,pointsGFMC, markersize = sizes,colormap = :viridis, color = sizefunc.(corrEnd),alpha = 1.0,marker = markerfunc.(corrEnd))
    fig
end


function getFT(vals,allPlaqs,refPlaq)
    ri = [ SW.SVector(r .- refPlaq) for r in allPlaqs]

    k = trueMomenta(0,2pi,size(S,1))

    FTgen = [FTPlaq(ri,vals,SW.SA[kx,ky]) for (kx,ky) in Iterators.product(k,k)]
end
BBMat = getFT.(BBCorrelator,Ref(allPlaqs),Ref(refPlaq))
##
with_theme(theme_PiTicks()) do
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    # FT = PlaqSumFT(S,BBCorrelator,refPlaq,allPlaqs)

    # hm = heatmap!(ax,kx,ky,FT,colormap = :viridis)
    hm = heatmap!(ax,kx,ky,-mean(BBMat),colormap = :viridis)
    # hm = heatmap!(ax,kx,ky,BBMat[10],colormap = :viridis)
    Colorbar(fig[1,2],hm)
    fig
end

##
function getUncorrPart(S,refPlaq,k = trueMomenta(0,2pi,size(S,1))) 
    AllPlaqs = collect(SW.plaquetteIterator(S))

    ri = [ SW.SVector(r .- refPlaq) for r in AllPlaqs]

    FTGen = zeros(length(k),length(k))

    for (i,kx) in enumerate(k)
        for (j,ky) in enumerate(k)
            FT = 0.
            for r in ri
               FT += cos(kx*r[1] + ky*r[2])
            end
            FTGen[i,j] = FT
        end
    end
    FTGen/length(AllPlaqs)
end
##

S = SW.stencilConfig(zeros(12,12),1;boundaryCondition = :periodic)
S .= SW.periodicStateDenseLoops(size(S,1))
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.5))
CTSR = SW.ContinuousTimeMethod(10*CT.τ,Hxx = CT.Hxx)

ψG = SW.SimpleJastrowFunction(S)
# ψG = SW.RKFunction()
##
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,10,ψG,300,1e-3,SW.IterativeSRSolver();Nwalkers = 1*28,reconfigure = false,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=10_000)
SW.get_params(ψG) .= stochReconfRes.params
plotVarEn(stochReconfRes,movavg = 30)
##

qvals = let 
    qx = trueMomenta(0,pi,size(S,1))
    qy = trueMomenta(0,pi,size(S,2))
    qvals = [SA[qx,qy] for (qx,qy) in Iterators.product(qx,qy)][:]
    filter!(x->x[1]>=x[2],qvals)
end 
##
BOp = SW.RandomPlaquetteFlipOperator(S)
Bi = fetch.([Threads.@spawn SW.measureObservables(S,BOp,[nothing],20,20*2,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:15]) ./ length(collect(SW.plaquetteIterator(S)))
# resB = fetch.([Threads.@spawn SW.measure_operator(S,CT,res.SaveConfigs,1,BOp,ψG,collect(SW.plaquetteIterator(S))[1:1]) for (i,res) in enumerate(results)])
##
BBQOp = SW.BBqOperator()
BBQ_0 = fetch.([Threads.@spawn SW.measureObservables(S,BBQOp,[SA[0,0.]],20,20*2,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:15])

BBQ = fetch.([Threads.@spawn SW.measureObservables(S,BBQOp,qvals,20,20*2,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:15])

##
Bq = fetch.([Threads.@spawn SW.measureObservables(S,SW.BqOperator(),qvals,20,20*2,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:15])
##
function processBBQ(S,qvals,BBQ,Bi,BBq_0,p=lastindex(Bi))
    allPlaqs = collect(SW.plaquetteIterator(S))
    FTmean = BBQ[p,:] ./length(allPlaqs)
    # kx,ky,FTrec = makeMatrix(qvals,FTmean)
    kx,ky,FTrec = makeMatrix(reconstruct_momentumSpace(qvals,FTmean)...)
    # return FTrec
    FTrec .-= FTrec[1,1]*0.5
    # FTrec .-= BBq_0[p,1]*0.5./length(allPlaqs)
    Bimean = Bi[p,1]

    refPlaq = (1,2)
    Bq_sq = Bimean^2 .* getUncorrPart(S,refPlaq,kx) *length(allPlaqs)
    # Bq_sq = getFT(Bi[p,:].^2,allPlaqs,refPlaq)

    FTrec .-= Bq_sq
    # FTrec = Bq_sq
    # return (;kx,ky,BBq = -Bq_sq)
    return (;kx,ky,BBq = FTrec)
end
##
function getBq_square(qvals,Bq)
    kx,ky,BqRec = makeMatrix(reconstruct_momentumSpace(qvals,Bq)...)
    BqRec .-= BqRec[1,1]*0.5

    Bqfunc = SW.getSqCont(BqRec)
    Bqsquare = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            Bq_re = Bqfunc(kx,ky)^2
            Bq_im = Bqfunc(kx-pi,ky-pi)^2
            Bqsquare[i,j] = Bq_re + Bq_im
        end
    end
    return Bqsquare
end
function processBBQ2(S,qvals,BBQ,Bq,BBq_0,p=lastindex(Bi))
    allPlaqs = collect(SW.plaquetteIterator(S))
    FTmean = BBQ[p,:] ./length(allPlaqs)
    # kx,ky,FTrec = makeMatrix(qvals,FTmean)
    kx,ky,FTrec = makeMatrix(reconstruct_momentumSpace(qvals,FTmean)...)
    # kx,ky,FTrec = makeMatrix(qvals,FTmean)
    
    Bq_sq = getBq_square(qvals,Bq[p,:]) /length(allPlaqs)

    # FTrec .-= FTrec[1,1]*0.5
    FTrec .-= BBq_0[p,1]*0.5./length(allPlaqs)
    
    FTrec .-= Bq_sq

    return (;kx,ky,BBq = Bq_sq)
    return (;kx,ky,BBq = FTrec)
end

##
with_theme(theme_PiTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1];aspect = 1)
    resBq = processBBQ2.(Ref(S),Ref(qvals),BBQ,Bq,BBQ_0,10)
    # resBq = processBBQ.(Ref(S),Ref(qvals),BBQ,Bi,BBQ_0,10)
    # resBq = processBBQ.(Ref(S),Ref(qvals),BBQ,BVals,BBQ_0,10)
    kx,ky = resBq[1].kx,resBq[1].ky
    BBq = [BBq.BBq for BBq in resBq]
    return BBq
    # println(maximum(mean(BBq)))
    hm = heatmap!(ax,kx,ky,mean(BBq),colormap = :viridis)
    # hm = heatmap!(ax,mean(Bq),colormap = :viridis)
    # hm = heatmap!(ax,kx,ky,Bq_sq,colormap = :viridis)
    Colorbar(fig[1,2],hm)
    fig
end
##
##
function ω_photon(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
ω_photon((kx,ky)) = ω_photon(kx,ky)
