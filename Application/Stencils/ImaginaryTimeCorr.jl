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

S = SW.stencilConfig(parent(SW.getStairCase(8)),1/2;
boundary = SW.Stencils.Wrap(),padding = SW.Stencils.Conditional()
)
S_ED = getPeriodic(SW.getStairCase(size(S,1)))
# S_ED = SW.getStairCase(size(S,1))
# S = SW.stencilConfig(parent(SW.getStairCase(7)),1/2)
# S_ED = SW.getStairCase(size(S,1))
HStair = SW.generateHilbertSpace(S_ED)

ExSol = SW.eigen(Matrix(HStair.H))
E0 = ExSol.values[1]
v0 = ExSol.vectors[:,1]
HConfs = SW.spinConfig.(HStair.AllStates,Ref(S_ED),Ref(HStair.plaqMapping))
##
function constructExactGuidingFunc(v0,AllStates)
    AllSTDict = Dict(SW.stencilConfig(parent(s),1/2)=>i for (i,s) in enumerate(AllStates))
    function psiG(Conf::SW.StencilSpinConfig)
        ind = AllSTDict[Conf]
        return v0[ind] #+ 0.01
    end
end

##

nThermal = 1000
nBra = 3
# ψG = SW.PlaquetteNumberGuidingFunction(0.197)
ψG = constructExactGuidingFunc(v0,HConfs)
##
# ψG(N) = 1
CT = SW.ContinuousTimeMethod(0.1,1,-E0)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,40,2500,ψG,equilibration_steps=100,pre_equilibration_steps=1_000,scatter_fraction=0.5) for i in 1:32])
##
# plotEnergies(results,CT,E0,nThermal=10,dense=true)
plotEnergies(results,CT,E0,nThermal=10,dense = false,Emin = E0-2e-1,Emax = E0+1e-1,τ = 1)
# plotEnergies(results,CT,E0,nThermal=10,normalize=true,Emin=E0-2e-2,Emax = E0+2e-2,dense = true)
##
function getSStau(res,p)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    ObsFunc = SW.constructSpinCorrTau(res.SaveConfigs)

    SW.getImagTimeCorr(Gnp,res.reconfigurationTable,ObsFunc)
end

AllCorrs = getSStau.(results,100)
meanSStau = mean(AllCorrs)
##
tau = CT.τ .* (eachindex(meanSStau) .-1)
S_tau_exact = SW.getSiSi_tau(HConfs,ExSol,tau)
##
lines(tau,meanSStau,label = "Spin Correlation",color = :black)
lines!(tau,S_tau_exact,label = "Spin Correlation",color = :red,linestyle = :dash)
ylims!(0,1/4)
band!(tau,meanSStau .- sqrt.(std(AllCorrs)),meanSStau .+ sqrt.(std(AllCorrs)),color = (:black,0.3))
current_figure()