using SpiderWebModel
import SpiderWebModel as SW
using SpiderWebModel.JuMP
using SpiderWebModel.HiGHS
using MakieHelpers
using CairoMakie
using SpiderWebModel.Gurobi
using ProgressMeter
using CairoMakie.Makie.ColorSchemes
using HDF5
include("SpiderwebDefs.jl")
const GRB_ENV = Gurobi.Env()
##
function setGurobiParameters!(model;kwargs...)
    JuMP.set_attribute(model, "GomoryPasses", 0)
    # Heuristics: expensive but good to find feasible solutions
    JuMP.set_optimizer_attribute(model, "Heuristics", 0.5)
    JuMP.set_optimizer_attribute(model, "ZeroObjNodes", 100)
    # JuMP.set_optimizer_attribute(model, "PumpPasses", -1)
    # JuMP.set_optimizer_attribute(model, "MinRelNodes", 1000)
    # JuMP.set_attribute(model, "Method", 2)
    JuMP.set_attribute(model, "PrePasses", 2)
    JuMP.set_attribute(model, "MIPFocus", 1)
    JuMP.set_optimizer_attribute(model, "PoolSearchMode", 0)   # 1: non-exhaustive search mode
    JuMP.set_optimizer_attribute(model, "PoolSolutions", 1)
    for (k,v) in kwargs
        JuMP.set_optimizer_attribute(model, string(k), v)
    end

end

function setUpSpiderWeb(L;OutputFlag=false,TimeLimit=3*60,LogFile="ConfigSearch_$L.log")
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    setUpSpiderWeb!(model,L)

    if OutputFlag
        JuMP.set_optimizer_attribute(model, "LogFile", LogFile)
        JuMP.set_optimizer_attribute(model, "LogToConsole", 0)  # number of solutions
    end
    setGurobiParameters!(model;TimeLimit,OutputFlag)
    return model
end

function demoInitializer(L,initializeDenominator=12)
    Conf = SW.SpinConfig(fill(NaN,L,L),1/2)
    numFixed = L^2÷initializeDenominator
    fixInds = rand(CartesianIndices(Conf),numFixed)
    vals = rand(-0.5:0.5,numFixed)
    for (I,v) in zip(fixInds,vals)
        Conf[I] = v
    end
    SW.plotSpinConfig(Conf)
end

function constructGroundstate(model,L,numSols)
    Sz = model[:Sz]

    if numSols >1
        JuMP.set_optimizer_attribute(model, "PoolSearchMode", 1)   # 1: non-exhaustive search mode
    end

    JuMP.set_optimizer_attribute(model, "PoolSolutions", numSols)
    JuMP.optimize!(model)

    num_results = result_count(model)
    
    Results = BitMatrix[] 
    step_successful = is_solved_and_feasible(model)
    if !step_successful
        return Results,step_successful
    end
    for n in 1:num_results
        sol = BitMatrix(undef,L,L)
        for i in eachindex(Sz)
            sol[i] = round(Bool,JuMP.value(Sz[i]; result=n))
        end
        push!(Results,sol)
    end

    return Results,step_successful
end

function multipleRuns(L,modelfactory,setups,NIndiv)
    sols = BitMatrix[]
    p = Progress(length(setups))
    # x = [BitMatrix(undef,0,0)]
    
    generate_showvalues(iter, x) = () -> [(:iter,iter), (:numSolutions,length(x))]
    sols = Vector{Vector{BitMatrix}}(undef,length(setups))
    statuses = falses(length(setups))
    # for iter in eachindex(setups)
    Threads.@threads for iter in eachindex(setups)
        model = modelfactory(L)
        set_attribute(model, MOI.NumberOfThreads(), 1)
        # model = StartModel
        fixinds,vals = setups[iter]
        fixSpins!(model,fixinds,vals)

        x,status = constructGroundstate(model,L,NIndiv)
        # @info "" iter length(x) typeof(x)
        sols[iter] = x
        statuses[iter] = status
        next!(p; showvalues = generate_showvalues(iter, sum(statuses)))
    end
    return [s  for solsvec in sols for s in solsvec if !isempty(s)]
end
##
L = 70
modelfactory(L) = setUpSpiderWeb(L,OutputFlag=false,TimeLimit=20,LogFile="ConfigSearch_$L.log")
setups = [getRandomSpins(L,8) for _ in 1:1000]
# JuMP.optimize!(model)
##
sols = multipleRuns(L,modelfactory,setups,1)
##
oldsols = h5read("../ConfsRaw/Configs_40.h5", "Confs")
Confs = SW.floatSpinConfig.(SW.SpinConfig.(eachslice(oldsols,dims=3),1/2))
# newsols = cat(oldsols,Array{Bool,3}(stack(sols,dims=3)),dims=3)
##  1
# rm("Configs_$L.h5")
# h5write("Configs_$L.h5", "Confs", newsols)
##
Confs = SW.floatSpinConfig.(SW.SpinConfig.(sols,1/2))
##
allsols = cat(h5read("../ConfsRaw/Configs_40.h5","Confs"),h5read("../ConfsRaw/Configs_40.h5","Confs2_6"),dims=3)
# allsols = h5read("../ConfsRaw/Configs_40.h5","Confs")
Confs = SW.floatSpinConfig.(SW.SpinConfig.(eachslice(allsols,dims=3),1/2))
##
SW.plotFractons(rand(Confs))
##
# SqLargeN(qx,qy) = -((2*(cos(qx) - cos(qy) + 2sin(qx)*sin(qy))^2)/(-4 + 
SqLargeN(qx,qy) = 2 - (cos(qx) - cos(qy))^2/((cos(qx) - cos(qy))^2 + (-cos(qx - qy) + cos(qx + qy))^2) - (2*(cos(qx) - cos(qy))*(-cos(qx - qy) + cos(qx + qy)))/((cos(qx) - cos(qy))^2 + (-cos(qx - qy) + cos(qx + qy))^2) - (-cos(qx - qy) + cos(qx + qy))^2/((cos(qx) - cos(qy))^2 + (-cos(qx - qy) + cos(qx + qy))^2)

Sq = SW.getEqualWeightStructureFac(Confs)
with_theme(theme_PiTicks()) do
    Sqre = real(Sq.Sq_k)

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)
    L = size(Confs[1],1)
    ferroVal = real(Sq.Sq_k[1,1])
    
    maxVal = argmax(≠(ferroVal), Sqre)
    @info "" sumrule = sum(real(Sq.Sq_k[1:end-1,1:end-1])) sumRuleExpected = L^2/4

    # ratio_viridis = maxVal / ferroVal +0.3
    # ratio_greys3 = 1 - ratio_viridis
    # uniquevals = unique(Sqre)

    # viridis = get(ColorSchemes.viridis, LinRange(0,maxVal,round(Int,ratio_viridis*100)))
    # greys3 = get(cgrad([:grey,:red]), LinRange(0,ferroVal,round(Int,ratio_greys3*100)))
    # virgre3 = vcat(viridis, greys3)
    # vals = collect(LinRange(0,ferroVal^(1.9),lenCol).^(1/1.9))
    # push!(vals,ferroVal)
    # hm = heatmap!(ax,Sq.k,Sq.k,Sqre)
    SqLN = [qx > pi ? SqLargeN(qx,qy)/2 : NaN for qx in Sq.k, qy in Sq.k]
    hm = heatmap!(ax,Sq.k,Sq.k,Sqre)
    # heatmap!(ax,Sq.k,Sq.k,SqLN,colorrange=(0,1.5))
    # heatmap!(ax,Sq.k,Sq.k,[qx == qy ? 1 : NaN for qx in Sq.k, qy in Sq.k],colormap=cgrad([:white]))
    # lines!(ax,Sq.k,Sq.k,color = :white,lw = 1)
    # lines!(ax,[pi-1e-10,pi+1e-10],[-0.2,2.2pi],color = :white,lw = 1)
    # hm = heatmap!(ax,Sq.k,Sq.k,Sqre,colorrange=(0,ferroVal),colormap=:viridis)
    
    # hms = [S == ferroVal ? ferroVal : NaN for S in Sqre]
    # heatmap!(ax,Sq.k,Sq.k,hms,colorrange=(0,ferroVal),)
    # push!(colors,ColorSchemes.RGB(1.,0,0))
    xlims!(ax,0,2pi)
    ylims!(ax,0,2pi)
    Colorbar(fig[1,2],hm)
    fig

    # save("Application/figs/Sq_$L.pdf",fig)
end
