using JuMP
# import ConstraintSolver as CS
import HiGHS
include("../SpiderwebDefs.jl")
using Random
import Gurobi
HiGHS.Highs_resetGlobalScheduler(1)

##
Random.seed!(12345)
L = 40
setups = [getRandomSpins(L,6) for i in 1:10]
##
function runsetups(model,setups)
    times = Float64[]
    solved = Bool[]
    solutions = BitMatrix[]
    for (fixInds,vals) in setups
        fixSpins!(model,fixInds,vals)
        t = @elapsed JuMP.optimize!(model)
        status = is_solved_and_feasible(model)
        if status
            sol = BitMatrix(round.(Bool,JuMP.value.(model[:y])))
            push!(solutions,sol)
        end
        push!(solved,status)
        push!(times,t)
    end
    return (;times,solved,solutions)
end

function setAttributes(model,solver::Type{Gurobi.Optimizer})
    set_attribute(model, MOI.NumberOfThreads(), 1)
    #setting presolver flags to be aggressive 
    # JuMP.set_optimizer_attribute(model, "Presolve", 2) #aggressive
    # JuMP.set_optimizer_attribute(model, "PrePasses", 2_000_000_000) #many prepasses
    # JuMP.set_optimizer_attribute(model, "Aggregate", 2) #aggressive
    # JuMP.set_optimizer_attribute(model, "AggFill", 2_000_000_000) #aggressive
    # JuMP.set_optimizer_attribute(model, "PreSparsify", 2) #aggressive
    
    JuMP.set_attribute(model, "OutputFlag", false)
    JuMP.set_attribute(model, "GomoryPasses", 0)
    JuMP.set_optimizer_attribute(model, "PoolSolutions", 1)
    
    # Heuristics: expensive but good to find feasible solutions
    JuMP.set_optimizer_attribute(model, "Heuristics", 0.5)
    JuMP.set_optimizer_attribute(model, "ZeroObjNodes", 100)
    # JuMP.set_optimizer_attribute(model, "PumpPasses", -1)
    # JuMP.set_optimizer_attribute(model, "MinRelNodes", 1000)

    # JuMP.set_attribute(model, "Method", 2)
    JuMP.set_attribute(model, "PrePasses", 2)
    JuMP.set_attribute(model, "MIPFocus", 1)
    JuMP.set_optimizer_attribute(model, "IntegralityFocus", 1) #prioritize integrality, i.e. integer solutions
    # grb = backend(model)
    # JuMP.set_optimizer_attribute(solver, "TimeLimit", 30)
end

function setAttributes(model,solver::Type{HiGHS.Optimizer})
    set_attribute(model, MOI.NumberOfThreads(), 1)
    # JuMP.set_attribute(model, "OutputFlag", false)
end

setAttributes(_,_) = nothing # default fallback

function runSolver(solver,setups)
    model = Model(solver)
    set_time_limit_sec(model, 30.0)
    setUpSpiderWeb!(model,L)
    setAttributes(model,solver)
    runsetups(model,setups)
end

##
# statsHighs = runSolver(HiGHS.Optimizer)
##
@time statsGur = runSolver(Gurobi.Optimizer,setups[1:1])
println("num solutions:" , sum(statsGur.solved))
