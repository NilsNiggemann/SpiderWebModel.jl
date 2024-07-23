import Pkg
cd(@__DIR__)
Pkg.activate("../.")
import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel
using Optim
include("plottingUtils.jl")
##
#___________Periodic Boundaries_______________________
function getPeriodic(parent)
    state = parent |> Array
    SW.SpinConfig(SW.PeriodicMatrix(state), parent.S)
end

S = SW.stencilConfig(parent(SW.getStairCase(10)),1/2;
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
ψG = SW.PlaquetteNumberGuidingFunction(0.197)
# ψG = constructExactGuidingFunc(v0,HConfs)
##
# ψG(N) = 1
CT = SW.ContinuousTimeMethod(0.025,1,-E0)
# @time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,2000,52000,ψG,equilibration_steps=500,pre_equilibration_steps=1_000,scatter_fraction=0.5,
# outfile = "../../temp/$(i)_4.h5"
# ) for i in 1:24])
@time println("done")
##
results = vcat(SW.readResults.(filter(contains("_4.h5"),readdir("../../temp",join=true)),32000)...)
# plotEnergies(results,CT,E0,nThermal=10,dense=true)
plotEnergies(results,CT,E0,nThermal=10,dense = true,Emin = E0-2e-3,Emax = E0+1e-3,τ = 50)
# plotEnergies(results,CT,E0,nThermal=10,normalize=true,Emin=E0-2e-2,Emax = E0+2e-2,dense = true)
##
function getSStau(res,p,mtau,O)
    Gnp = SW.precomputeNormalizedAccWeight(res.TotalWeights,1,p)
    ObsFunc = SW.constructSWF_operator(res.SaveConfigs,O)
    
    GSObs = SW.getObs(Gnp,res.SaveConfigs,res.reconfigurationTable,O)

    # O = SW.PlaquetteNumberOperator_1(I)
    # function ObsFunc(α,n)
    #     conf = @view res.SaveConfigs[:,:,α,n]
        
    #     confS = SW.StencilSpinConfig(
    #         SW.Stencils.StencilArray{Tuple{8,8}}(conf,SW.Stencils.Moore(),SW.Stencils.Wrap(),SW.Stencils.Conditional()),
    #         Int8(1)
    #     )
        
    #     return O(confS)

    # end

    corrs = SW.getImagTimeCorr(Gnp,res.reconfigurationTable,ObsFunc,mtau) 
    for i in eachindex(corrs)
        corrs[i] .-= GSObs.^2
    end
    meancorrs = mean.(corrs) 
end
plaqs = collect(SW.plaquetteIterator(S))

# AllCorrs = [getSStau(res,100,I) for I in plaqs for res in results]
# @profview getSStau.(results[1:4],10,Ref((2,3)))
# obs(x) = x[2,3]+x[3,3]#*x[2,3]
obs(x,normalization=1) = x[2,3]*normalization^2*x[2,4] + x[2,3]*normalization
function ConstructObsMat(S,normalization=1.)
    # normConv = convert(Float32,normalization)
    # xFl = zeros(Float32,size(S))
    # # x -> x .*normalization
    # x -> xFl .= x.*normConv
    # identity
    c = copy(S)
    function obs(x)
        c .= x
        SW.NPlaquettes(x)
    end
end
#  AllCorrs = getSStau.(results,20,x->obs(x),0.5)
# @profview getSStau.(results,20,x->obs(x,0.5))
# AllCorrs = fetch.([Threads.@spawn getSStau(res,200,30,ConstructObsMat(S,0.5)) for res in results])
@time AllCorrs = fetch.([Threads.@spawn getSStau(res,650,160,identity) for res in results])
# AllCorrs = fetch.([Threads.@spawn getSStau(res,200,30,x->obs(x,0.5)) for res in results])
# AllCorrs = getSStau.(results,200,x->x[1,2]*0.5 + x[2,3]*0.5)
meanSStau = mean(AllCorrs) .*0.25 #1/4 since we compute <S_i(tau) S_i(0)> and spins are rescaled by a factor of 2
##
tau = CT.τ .* (eachindex(meanSStau) .-1)
# S_tau_exact = SW.getTauCorr(HConfs,ExSol,tau,obs) .- SW.getGSObsED(HConfs,v0,obs)^2

S_tau_exact =  let
    S_tau_exact_Mat = SW.getTauCorr(HConfs,ExSol,tau,identity)
    a = SW.getGSObsED(HConfs,v0,identity).^2
    for i in eachindex(S_tau_exact_Mat)
        S_tau_exact_Mat[i] .-= a
    end
    mean.(S_tau_exact_Mat)
end
##
function fitExp(tau,O_tau,numTerms=1;verbose=false)

    function model(A_n,En,tau)
        return sum(a^2 * exp(-abs(e) * tau) for (a,e) in zip(A_n,En))
    end
    model(v,tau) = @views model(v[1,:],v[2,:],tau)

    function loss(params)
        predic = model.(Ref(params),tau)
        return sum(abs2,(O_tau .- predic) ./O_tau)
    end


    x0 = ones(2,numTerms)
    # res = optimize(loss, x0,LBFGS(); autodiff = :forward)
    # precond(n) = SW.SparseArrays.spdiagm(-1 => -ones(n-1), 0 => 2ones(n), 1 => -ones(n-1)) * (n+1)

    # res = optimize(loss, x0, method = ConjugateGradient(P = nothing))
    res = optimize(loss, x0)
    # res = optimize(loss, x0,SimulatedAnnealing())
    
    verbose && display(res)
    sol = Optim.minimizer(res)
    A_i = sol[1,:]
    Δ_i = sol[2,:]
    return (;A_i,Δ_i,sol,model,res)
end
##
with_theme(theme_SimpleTicks()) do 
    nfit = 4
        
    firstInd = findfirst(>(0.),tau)
    lastInd = findfirst(>=(4.8),tau)

    fig = Figure(size = 0.7 .*(800, 1000))
    ax = Axis(fig[1, 1], xlabel = L"τ", ylabel = L"\mathcal{O}(τ)",
    # yscale = Makie.pseudolog10,
    yscale = log10,
    xlabelvisible = false, xticklabelsvisible=false
    )

    axfit = Axis(fig[2, 1], xlabel = L"τ", ylabel = L"|\mathcal{O}(τ) - \mathcal{O}_\textrm{fit}(τ)|",yscale = log10,xlabelvisible = false, xticklabelsvisible=false)

    axdelta = Axis(fig[3, 1], xlabel = L"τ", ylabel = L"\Delta")

    if isnothing(lastInd)
        lastInd = lastindex(tau)
    end
    if isnothing(firstInd)
        firstInd = firstindex(tau)
    end
    # (;A,Delta)= fitExp(tau[firstInd:end],meanSStau[firstInd:end])
    taufit = tau[firstInd:lastInd]
    meanSStaufit = meanSStau[firstInd:lastInd]

    (;A_i,Δ_i,sol,model) = fitExp(taufit,meanSStaufit,nfit)
    
    println(sort(Δ_i))
    
    Delta = minimum(Δ_i)
    cutoff = 1e-10

    normalizeVals(x) = max(x,cutoff)
    yGFMC = normalizeVals.(meanSStau)
    yExact = normalizeVals.(S_tau_exact)
    yGFMCUpper = normalizeVals.(meanSStau .+ std(AllCorrs))
    yGFMCLower = normalizeVals.(meanSStau .- std(AllCorrs))

    scatterlines!(ax,tau,yGFMC,label = L"GFMC $$",color = :black)
    lastval = S_tau_exact[end] 
    lastval_tau = taufit[end]

    gap = (ExSol.values[2] - ExSol.values[1])

    y0 = lastval * exp(gap * lastval_tau)
    lines!(ax,taufit,normalizeVals.(model.(Ref(sol),taufit)),label = L"fit $\sum_n^{%$nfit} |A_n|^2 e^{-\Delta_n \tau}$",color = :blue,linestyle = :solid,linewidth = 3)

    lines!(ax,tau,yExact,label = L"exact $$",color = :red,linestyle = :dash)

    
    # lines!(ax,taufit,y0 .*exp.(-taufit*gap),label = L"e^{ΔE}",color = :blue,linestyle = :dash)

    # lines!(ax,taufit,A .*exp.(-taufit*Delta),label = L"e^{ΔE}",color = :blue,linestyle = :dash)
    # ylims!(0,1/4)
    band!(ax,tau,yGFMCLower,yGFMCUpper,color = (:black,0.3))

    lowerlim = minimum(filter(x->x >cutoff+1e-6, yGFMC))

    ylims!(ax,lowerlim,maximum(yGFMC))

    Δs = [log.(normalizeVals.(O_GFMC[1])./normalizeVals.(O_GFMC)) ./ tau for O_GFMC in AllCorrs]
    Δ = mean(Δs)
    Δerr = std(Δs)
    scatterlines!(axfit,taufit,abs.(meanSStaufit .- model.(Ref(sol),taufit)),label = L"\mathcal{O}(τ) - \mathcal{O}_\textrm{fit}(τ)",color = :black)

    scatterlines!(axdelta,tau[firstInd:lastInd],Δ[firstInd:lastInd],label = L"\Delta_\textrm{GFMC}(τ)",color = :black)
    upperband = Δ .+ Δerr
    lowerband = Δ .- Δerr
    band!(axdelta,tau[firstInd:lastInd],upperband[firstInd:lastInd],lowerband[firstInd:lastInd],color = (:black,0.3))
    
    strd(x) = string(round(x,digits = 3))
    Δex = log.(yExact[1]./yExact) ./ tau
    lines!(axdelta,tau,Δex,label = L"\Delta_\textrm{exact}(τ)",color = :red,linestyle = :dash)
    hlines!(axdelta,[gap],color = :red,label = L"Δ_\textrm{exact} = %$(strd(gap))")
    hlines!(axdelta,[Delta],color = :blue,label = L"Δ_\textrm{fit} = %$(strd(Delta))")
    axislegend(ax)
    axislegend(axdelta,position=:rt)
    
    rowsize!(fig.layout, 1, Relative(0.4))
    rowsize!(fig.layout, 2, Relative(0.2))
    rowsize!(fig.layout, 3, Relative(0.4))
    linkxaxes!(ax, axfit, axdelta)
    save("Otau4.pdf", fig)
    fig
end

##
obs(x,normalization=1) = x[1,2]+ x[2,3] + x[4,5]#*normalization^3*x[5,4]#*x[1,3] #+ x[2,3]*normalization

with_theme(theme_SimpleTicks()) do 
    tau = LinRange(0,20.2,100)
    S_tau_exact = SW.getTauCorr(HConfs,ExSol,tau,obs) .- SW.getGSObsED(HConfs,v0,obs)^2

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"τ", ylabel = L"\mathcal{O}(τ)",
    yscale = log10,
    )
    axdelta = Axis(fig[2, 1], xlabel = L"τ", ylabel = L"\Delta")

    lines!(ax,tau,S_tau_exact,label = L"\mathcal{O}(τ)",color = :red)

    Δex = log.(S_tau_exact[1]./S_tau_exact) ./ tau
    lines!(axdelta,tau,Δex,label = L"\Delta_\textrm{exact}(τ)",color = :red,linestyle = :dash)
    gap = (ExSol.values[2] - ExSol.values[1])
    strd(x) = string(round(x,digits = 3))
    hlines!(axdelta,[gap],color = :red,label = L"Δ_\textrm{exact} = %$(strd(gap))")
    fig
end