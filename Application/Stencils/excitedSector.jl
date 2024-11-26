import Pkg

import SpiderWebModel as SW
using CairoMakie
using Statistics
using MakieHelpers
using SpiderWebModel

include("plottingUtils.jl")
include("4x4Orders_base.jl")
##
function shuffleSector!(S,N;maxiter = 10N)

    flips = (:diag,:anti)
    # flips = (:row,:col,:diag,:anti)
    spin = SW.getSpin(S)
    function transformSpins!(vec,sgn)
        any(==(sgn*spin),vec) && return
        vec .+= sgn
        return
    end

    function applyMove!(S)
        whichflip = rand(flips)
        sgn = -2sign(sum(S))
        if sgn == 0 
            sgn = rand((-2,2))
        end
        if whichflip == :row
            i = rand(axes(S,1))
            row = @view S[i,begin+iseven(i):2:end]
            transformSpins!(row,sgn)
        elseif whichflip == :col
            j = rand(axes(S,2))
            col = @view S[begin+iseven(j):2:end,j]
            transformSpins!(col,sgn)
        elseif whichflip == :diag
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,1,true)
            transformSpins!(diag,sgn)
        elseif whichflip == :anti
            i = rand(axes(S,2)[1:2:end])
            diag = SW.getDiagonal(S,i,-1,true)
            transformSpins!(diag,sgn)
        end
        return 
    end

    for iteration in 1:N
        applyMove!(S)
    end
    for iter in N:maxiter 
        applyMove!(S)
        sum(S) == 0 && break
    end
    # S .= 2S
    SW.fulFillsConstraint(S) || @warn "Constraint not fulfilled"
    return S
end
##
SW.Random.seed!(7870)
S = SW.stencilConfig(zeros(20,20),1;
boundaryCondition = :periodic
)
S .= SW.get4x4PeriodicState(size(S,1),20)
# SW.plotApplPlaquettes(S)

# SW.Random.seed!(1314)

shuffleSector!(S,5585400)
SW.plotApplPlaquettes(S)
##
ψG = SW.SimpleJastrowFunction(S,Float64)
SW.rand!(SW.getNonSymmetric(ψG),1e-3)
ψG.v_ij .= SW.Symmetric(ψG.v_ij)
CT = SW.ContinuousTimeMethod(0.1,Hxx = SW.Hxx_RK(0.7))
CTSR = SW.ContinuousTimeMethod(10,Hxx = CT.Hxx)
##
function stepSize(i)
    i > 500 && return 8e-4
    i > 200 && return 3e-4
    i > 100 && return 1e-4
    i > 1 && return 5e-5
    return 3e-3
end

SW.Random.seed!(134)
stochReconfRes = SW.stochastic_reconfiguration(S,CTSR,30 ,ψG,1000,2e-4,SW.IterativeSRSolver();Nwalkers = 1*20,reconfigure=false,rel_tolerance=0,equilibration_steps=1000,pre_equilibration_steps=40_000,
report_steps = 5,
reset = false,
)

ψGnew = deepcopy(ψG)

SW.get_params(ψGnew) .= stochReconfRes.params

plotVarEn(stochReconfRes,movavg = 20,alpha_index = 1)

##
nThermal = 1000
NWalkers = 10*20
NSteps = 5000
SW.Random.seed!(1234)
@time results = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,ψGnew;equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:20])
# @time results2 = fetch.([Threads.@spawn SW.startManyWalkerGFMC(S,CT,NWalkers,NSteps,SW.PlaquetteNumberGuidingFunction(0.2);equilibration_steps=nThermal,pre_equilibration_steps=nThermal) for _ in 1:20])
##
# plotEnergies(results2,CT;normalize=true,dense=true,τ = 10,color = :black)
plotEnergies(results,CT;normalize=true,dense=true,τ = 20,color = :red)
current_figure()
##
Sqs = SW.getSqsGFMC(results,1:150)

##
with_theme(theme_PiTicks()) do 
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1, xlabel = L"q_x",ylabel = L"q_y")
    k = trueMomenta(0,2pi,size(S,1))
    hm = heatmap!(ax,k,k,dropmean(Sqs,dims=4)[:,:,120],colormap = :viridis)
    Colorbar(fig[1,2],hm)
    fig
    
end