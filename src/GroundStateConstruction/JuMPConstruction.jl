function setUpSpiderWeb(optimizer,L;)
    model = JuMP.Model(optimizer)
    setUpSpiderWeb!(model,L,optimizer())
end

function setUpSpiderWeb!(model,L,opt=nothing)
    INDEX = 1:L
    JuMP.@variable(model, Sz[INDEX,INDEX], Bin)
    setConstraints!(model,L,opt)
end

function setConstraints!(model,L,opt=nothing)
    Sz = model[:Sz]
    for i in 1:L
        for j in 1:L
            iseven(i+j) || continue
            plaquetteIsInBounds(Sz, i, j) || continue
            setConstraint!(model,Sz,i,j,opt)
        end
    end
    return model
end

function setConstraint!(model,Sz,i,j,opt)

    JuMP.@constraint(model,
    +Sz[i, j + 1]
    +Sz[i - 1, j + 1]
    -Sz[i - 1, j]
    -Sz[i - 1, j - 1]
    +Sz[i, j - 1]
    +Sz[i + 1, j - 1]
    -Sz[i + 1, j]
    -Sz[i + 1, j + 1] == 0)
end


function fixSpins!(model,fixInds,vals)
    
    Sz = model[:Sz]

    for i in eachindex(Sz)
        if JuMP.is_fixed(Sz[i])
            JuMP.unfix(Sz[i])
        end
    end
    for (I,v) in zip(fixInds,vals)
        JuMP.fix(Sz[I], v; force = true)
    end
    return model
end

function getRandomSpins(L::Integer,FixedFraction::Real)
    # numFixed = L^2 ÷ FixedFraction
    numFixed = round(Int,L^2*FixedFraction)
    fixInds = unique!(rand(CartesianIndices((L,L)),numFixed))
    vals = rand(0:1,length(fixInds))
    return fixInds,vals
end

function setUpSpiderWeb(L;kwargs...)
    model = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV_REF[]))
    JuMP.set_silent(model)
    setUpSpiderWeb!(model,L)
    # JuMP.set_optimizer_attribute(model, "LogFile", LogFile)
    setGurobiParameters!(model;OutputFlag=false,TimeLimit=3*60,kwargs...)
    return model
end


function setGurobiParameters!(model;kwargs...)
    JuMP.set_attribute(model, "GomoryPasses", 0)
    JuMP.set_optimizer_attribute(model, "LogToConsole", 0)  # number of solutions

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

function constructGroundstate(model,L)
    Sz = model[:Sz]

    JuMP.optimize!(model)

    step_successful = JuMP.is_solved_and_feasible(model)
    
    if !step_successful
        return BitMatrix(undef,0,0),step_successful
    end
    
    sol = BitMatrix(undef,L,L)
    for i in eachindex(Sz)
        sol[i] = round(Bool,JuMP.value(Sz[i]))
    end

    return sol,step_successful
end

function constructGroundstates(L::Integer,modelfactory::Function,NRuns::Integer,fixedFraction::Real)
    p = ProgressMeter.Progress(NRuns)
    
    generate_showvalues(iter, x) = () -> [(:iter,iter), (:numSolutions,length(x))]
    sols = Vector{BitMatrix}(undef,NRuns)
    statuses = falses(NRuns)
    Threads.@threads for iter in 1:NRuns
        model = modelfactory(L)
        JuMP.set_attribute(model, MOI.NumberOfThreads(), 1)

        fixinds,vals = getRandomSpins(L,fixedFraction)
        fixSpins!(model,fixinds,vals)
        
        x,status = constructGroundstate(model,L)
        sols[iter] = x
        statuses[iter] = status
        ProgressMeter.next!(p; showvalues = generate_showvalues(iter, sum(statuses)))
    end
    filter!(!isempty,sols)
end


function constructGroundstates(L::Integer,NRuns::Integer,fixedFraction::Real;kwargs...)
    modelfactory(L) = setUpSpiderWeb(L,OutputFlag=false,TimeLimit=20;kwargs...)
    constructGroundstates(L,modelfactory,NRuns,fixedFraction::Real)
end