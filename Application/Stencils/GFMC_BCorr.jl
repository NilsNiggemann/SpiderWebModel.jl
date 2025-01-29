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
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,40,ψG,100,3e-3,SW.IterativeSRSolver();Nwalkers = 20,rel_tolerance=1e-8,equilibration_steps=nThermal,pre_equilibration_steps=40_000,report_steps=10)
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
allPlaqs = collect(SW.plaquetteIterator(S))
active_Plaqs = SW.getApplicablePlaquettes(S_ED)
refPlaq = first(active_Plaqs)
# pairPlaqs = filter(!=(refPlaq),allPlaqs)
pairPlaqs = copy(active_Plaqs)
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
BVals = stack(fetch.([Threads.@spawn SW.measureObservables(S,BOp,active_Plaqs,5,20*1,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:10]))
##
BOp_rand = SW.RandomPlaquetteFlipOperator(S)
BValsRand = stack([SW.measureObservables(S,BOp_rand,[nothing],5,20*3,3000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:10]) ./ length(active_Plaqs)
##
# BVals = fetch.([Threads.@spawn SW.measureObservables(S,BOp,active_Plaqs,5,20*1,1000,CT,ψG;equilibration_steps = 3000,pre_equilibration_steps=20_000) for _ in 1:20])
##
BBOp = SW.BBOperator(S,refPlaq)
# GFMCPlaqs = SW.getApplicablePlaquettes(S)
GFMCPlaqs = active_Plaqs

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

function PlaqSumFT(S,BBCorr,refPlaq,allPlaqs,phase = 0)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                rx,ry = P .- refPlaq
                BBq += BBCorr[iP] * cos(kx*rx + ky*ry-phase)
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

function PlaqSumFT_complex(S,BBCorr,refPlaq,allPlaqs)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(ComplexF64,length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                rx,ry = P .- refPlaq
                BBq += BBCorr[iP] * exp(1im*(kx*rx + ky*ry))
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

function PlaqSumFTPairs(S,BBCorr::AbstractMatrix,allPlaqs,phase = 0)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                for (jP,P2) in enumerate(allPlaqs)
                    rx,ry = P .- P2
                    BBq += BBCorr[iP,jP] * cos(kx*rx + ky*ry-phase)
                end
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

function Cos_2_sum(S,BBCorr::AbstractMatrix,allPlaqs,phase = 0)
    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))
    FTres = zeros(length(kx),length(ky))
    for (i,kx) in enumerate(kx)
        for (j,ky) in enumerate(ky)
            BBq = 0.
            for (iP,P) in enumerate(allPlaqs)
                for (jP,P2) in enumerate(allPlaqs)
                    rx,ry = P .- P2

                    kr = kx*rx + ky*ry
                    BBq += BBCorr[iP,jP] * cos(0.5kr+phase)^2
                end
            end
            FTres[i,j] = BBq
        end
        
    end
    return FTres
end

## ____________BBQ Operator_______________________


BBQOp = SW.BBqOperator()
qvals = let 
    qx = trueMomenta(0,2pi,size(S,1))[1:end-1]
    qy = trueMomenta(0,2pi,size(S,2))[1:end-1]
    qvals = [SA[qx,qy] for (qx,qy) in Iterators.product(qx,qy)]
end
BBQOp = SW.BBqOperator()
##
BBQ = stack(fetch.([Threads.@spawn SW.measureObservables(S,BBQOp,qvals,5,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))

##
exactBB_pairs = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,Pi,Pj) for Pi in active_Plaqs, Pj in active_Plaqs] ./ length(allPlaqs)

with_theme(theme_PiTicks()) do 
    fig = Figure(size = (600,600))

    ax_exactB2q = Axis(fig[1,1];aspect = 1,title = L"exact $B^2_c(q)$")

    ax_MC_B2q = Axis(fig[2,1];aspect = 1,title = L"MC $B^2_c(q)$")

    ax_exactBBq = Axis(fig[1,3];aspect = 1,title = L"exact $\mathcal{B}(\mathbf{q})$")

    axMC_BBq = Axis(fig[2,3];aspect = 1,title = L"MC $\mathcal{B}(\mathbf{q})$")

    ax_diff_B2q_ex = Axis(fig[3,1];aspect = 1,title = L"diff $$")
    ax_diff_BBq_ex = Axis(fig[3,3];aspect = 1,title = L"diff $$")

    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    B2q_ex = Cos_2_sum(S,exactBB_pairs,active_Plaqs)

    BBq_ex = PlaqSumFTPairs(S,exactBB_pairs,active_Plaqs)
    # BBq_ex = 2B2q_ex .- B2q_ex[1,1]
    
    
    hm = heatmap!(ax_exactB2q,kx,ky,B2q_ex,colormap = :viridis)
    Colorbar(fig[1,2],hm)

    hm = heatmap!(ax_exactBBq,kx,ky,BBq_ex,colormap = :viridis)
    Colorbar(fig[1,4],hm)
    fig

    B2qmean = dropmean(BBQ,dims=3)[end,:]

    L_B2q= Int(sqrt(length(B2qmean)))

    B2q = SW.expand_Sq(reshape(B2qmean,L_B2q,L_B2q))

    hm = heatmap!(ax_MC_B2q,kx,ky,B2q,colormap = :viridis)
    Colorbar(fig[2,2],hm)
    
    hm = heatmap!(ax_diff_B2q_ex,kx,ky,B2q_ex-B2q,colormap = :viridis)
    Colorbar(fig[3,2],hm)
    BBq = 2B2q .-B2q[1,1]

    hm = heatmap!(axMC_BBq,kx,ky,BBq,colormap = :viridis)

    Colorbar(fig[2,4],hm)

    hm = heatmap!(ax_diff_BBq_ex,kx,ky,BBq_ex -BBq,colormap = :viridis)
    Colorbar(fig[3,4],hm)
    fig
end
##___________Bq Operator test_______________________

BQOp = SW.BqOperator()
BQ = stack(fetch.([Threads.@spawn SW.measureObservables(S,BQOp,qvals,5,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))

##
BQOp_sin = SW.BqSinOperator()
BQSin = stack(fetch.([Threads.@spawn SW.measureObservables(S,BQOp_sin,qvals,5,20*1,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:6]))
##
with_theme(theme_PiTicks()) do 
    fig = Figure(size = (600,600))
    axBqEx = Axis(fig[1,1];aspect = 1,title = L"\langle B^\textrm{re}(\mathbf{q})\rangle_\textrm{exact}")
    axBqSin = Axis(fig[1,3];aspect = 1,title = L"\langle B^\textrm{im}(\mathbf{q})\rangle_\textrm{exact}")

    axBqMC = Axis(fig[2,1];aspect = 1,title = L"\langle B^\textrm{re}(\mathbf{q})\rangle_\textrm{GFMC}")
    axBqSinMC = Axis(fig[2,3];aspect = 1,title = L"\langle B^\textrm{im}(\mathbf{q})\rangle_\textrm{GFMC}")

    axDiffBq = Axis(fig[3,1];aspect = 1,title = L"\Delta \langle B^\textrm{re}(\mathbf{q})\rangle")
    axDiffBqSin = Axis(fig[3,3];aspect = 1,title = L"\Delta \langle B^\textrm{im}(\mathbf{q})\rangle")

    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    # FTex = PlaqSumFT(S,exactB,SW.getCentralPlaquette(S),active_Plaqs,0) ./sqrt(length(allPlaqs))

    Bq_ex = PlaqSumFT_complex(S,exactB,SW.getCentralPlaquette(S),active_Plaqs) ./sqrt(length(allPlaqs))

    hm = heatmap!(axBqEx,kx,ky,real(Bq_ex),colormap = :viridis)
    Colorbar(fig[1,2],hm)

    # FTexSin = PlaqSumFT(S,exactB,SW.getCentralPlaquette(S),active_Plaqs,pi/2)./sqrt(length(allPlaqs))
    # FTexSin = PlaqSumFT_cos2(S,exactB,SW.getCentralPlaquette(S),active_Plaqs,pi/2)./sqrt(length(allPlaqs))

    hm = heatmap!(axBqSin,kx,ky,imag(Bq_ex),colormap = :viridis)
    Colorbar(fig[1,4],hm)

    
    Bqmean = dropmean(BQ,dims=3)[end,:]
    L_Bq= Int(sqrt(length(Bqmean)))
    B2q = SW.expand_Sq(reshape(Bqmean,L_Bq,L_Bq))

    # Bq = 2B2q .- B2q[1,1]
    Bq = SW.recover_Bq_real(B2q)

    hm = heatmap!(axBqMC,kx,ky,Bq,colormap = :viridis)
    Colorbar(fig[2,2],hm)

    BQSinmean = dropmean(BQSin,dims=3)[end,:]
    L_Bq= Int(sqrt(length(BQSinmean)))
    BqSin = SW.expand_Sq(reshape(BQSinmean,L_Bq,L_Bq))

    # BqSin = 2BqSin .- Bq[1,1]
    BqSin = SW.recover_Bq_imag(BqSin)
    # BqSin = 2BqSin .- 2BqSin[1,1]

    hm = heatmap!(axBqSinMC,kx,ky,BqSin,colormap = :viridis)
    Colorbar(fig[2,4],hm)

    diffBq = real(Bq_ex) - Bq
    hm = heatmap!(axDiffBq,kx,ky,diffBq,colormap = :viridis)
    Colorbar(fig[3,2],hm)

    diffBqSin = imag(Bq_ex) - BqSin
    hm = heatmap!(axDiffBqSin,kx,ky,diffBqSin,colormap = :viridis)
    Colorbar(fig[3,4],hm)

    fig
end
## _________together_______________________
BiBj = [Bi*Bj for Bi in exactB, Bj in exactB] ./length(allPlaqs)
# exactBB_pairs_decor = [SW.getBij_square(HStair.AllStates,HStair.plaqMapping,v0,Pi,Pj) for Pi in allPlaqs, Pj in allPlaqs] ./ length(allPlaqs)


with_theme(theme_PiTicks()) do 
    fig = Figure(size = (400,600))
    ax1 = Axis(fig[1,1];aspect = 1,title = L"\langle \mathcal{B}(\mathbf{q})\rangle_\textrm{exact}")
    
    BBQmean = dropmean(BBQ,dims=3)[end,:]
    
    L_Bq= Int(sqrt(length(BBQmean)))

    BQmean = dropmean(BQ,dims=3)[end,:]

    BQSinmean = dropmean(BQSin,dims=3)[end,:]
    
    B2q = SW.recover_Bq_real(BBQmean)

    CorrGFMC = SW.expand_Sq(reshape(B2q .- abs2.(SW.recoverBq(BQmean,BQSinmean)),L_Bq,L_Bq))
    
    ax2 = Axis(fig[2,1];aspect = 1,title = "GFMC Bq op")

    kx = trueMomenta(0,2pi,size(S,1))
    ky = trueMomenta(0,2pi,size(S,2))

    Corr_ex = PlaqSumFTPairs(S,exactBB_pairs .- BiBj ,active_Plaqs)

    hm = heatmap!(ax1,kx,ky,Corr_ex,colormap = :viridis)
    Colorbar(fig[1,2],hm)

    hm = heatmap!(ax2,kx,ky,CorrGFMC,colormap = :viridis)
    Colorbar(fig[2,2],hm)
    fig
end

##
#___________Spin-1_______________________
# S = SW.stencilConfig(zeros(16,16),1;boundaryCondition = :periodic)
# S .= SW.periodicStateDenseLoops(size(S,1))
S = SW.get4x4PeriodicSpinConf(16,2)
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.7))
CTSR = SW.ContinuousTimeMethod(20*CT.τ,Hxx = CT.Hxx)

ψG = SW.SimpleJastrowFunction(S)
ψGSymm = SW.symmetrize(ψG,SW.TranslationalSymmetry([2,2],[2,-2]),S)
SW.rand!(ψGSymm,1e-5)

# ψG = SW.RKFunction()
##

stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,100,ψGSymm,1000,8e-3,SW.IterativeSRSolver();Nwalkers = 1*20,rel_tolerance=1e-8,equilibration_steps=1000,pre_equilibration_steps=10_000,report_steps=5)
SW.get_params(ψG) .= stochReconfRes.params
plotVarEn(stochReconfRes,movavg = 30)
##
# scatter(Point.(allPlaqs))
# scatter!(Point.(reducedPlaqs),color = :red)
# current_figure()
allPlaqs = collect(SW.plaquetteIterator(S))

##

qvals = let 
    qx = trueMomenta(0,pi,size(S,1))[1:end-1]
    qy = trueMomenta(0,pi,size(S,2))[1:end-1]
    qvals = [SA[qx,qy] for (qx,qy) in Iterators.product(qx,qy)][:]
    filter!(x->x[1]>=x[2],qvals)
end 
##
BBQOp = SW.BBqOperator()

BBQ = fetch.([Threads.@spawn SW.measureObservables(S,BBQOp,qvals,5,20*2,1000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:10])

##
Bq_Cos2 = fetch.([Threads.@spawn SW.measureObservables(S,SW.BqOperator(),qvals,5,20*5,2000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:10])

Bq_Sin2 = fetch.([Threads.@spawn SW.measureObservables(S,SW.BqSinOperator(),qvals,5,20*5,2000,CT,ψG;equilibration_steps = 5000,pre_equilibration_steps=20_000) for _ in 1:10])
##
function getB2Q(qvals,BBQ,Bq_Cos2,Bq_Sin2,p=last(axes(BBQ,1)))

    @assert first(qvals) == SA[0,0]

    BBQ = SW.recover_Bq_real(BBQ[p,:])
    Bq_re = SW.recover_Bq_real(Bq_Cos2[p,:])
    Bq_im = SW.recover_Bq_imag(Bq_Sin2[p,:])

    B2Q = BBQ .- abs2.(Bq_re) .- abs2.(Bq_im)

    kx,ky,B2Q = makeMatrix(reconstruct_momentumSpace(qvals,B2Q)...)

    return (;kx,ky,BBq = B2Q)
end

getB2Q_p(qvals,BBQ,Bq_Cos2,Bq_Sin2) = [getB2Q(qvals,BBQ,Bq_Cos2,Bq_Sin2,p) for p in axes(BBQ,1)]
##

with_theme(theme_PiTicks()) do 
    fig = Figure()
    ax1 = Axis(fig[1,1];aspect = 1, title = L"\mathcal{B}(\mathbf{q})")
    ax2 = Axis(fig[1,2];aspect = 1, title = L"err $$")
    ax3 = Axis(fig[2,1:2];xticks = SimpleTicks(),yticks = SimpleTicks(),ylabel = L"\mathcal{B}(0,0)",xlabel = L"τ")
    resBq = getB2Q_p.(Ref(qvals),BBQ,Bq_Cos2,Bq_Sin2)
    kx,ky = resBq[1][1].kx,resBq[1][1].ky
    BBq = [BBq[5].BBq for BBq in resBq]

    BBqmean = mean(BBq)
    BBqstd = std(BBq)
    colorrange = extrema([extrema(BBqmean)...,extrema(BBqstd)...])

    # hm1 = heatmapSq!(ax1, mean(BBq),k_range=(-0.5pi,1.5pi), colormap = :viridis;colorrange)
    # hm2 = heatmapSq!(ax2, std(BBq),k_range=(-0.5pi,1.5pi), colormap = :viridis;colorrange)

    hm1 = heatmap!(ax1, kx, ky, mean(BBq), colormap = :viridis;colorrange)
    hm2 = heatmap!(ax2, kx, ky, std(BBq), colormap = :viridis;colorrange)

    Colorbar(fig[1,3], hm2)
    
    Is = [CartesianIndex(1,1), argmax(BBqmean),
        argmin(BBqmean),CartesianIndex(1,4)]

    for I in Is
        BBq_i = [[BBq[i].BBq[I] for i in 1:5] for BBq in resBq]
        
        errlines!(ax3, CT.τ.* (0:4), mean(BBq_i), std(BBq_i))
    end
    fig
end
##
function ω_photon(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
ω_photon((kx,ky)) = ω_photon(kx,ky)
