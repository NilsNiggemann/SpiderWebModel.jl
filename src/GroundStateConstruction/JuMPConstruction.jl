function setUpSpiderWeb(optimizer, L; kwargs...)
    model = JuMP.Model(optimizer)
    setUpSpiderWeb!(model, L, optimizer();kwargs...)
end

function setUpSpiderWeb!(model, L, opt = nothing;S=0.5, boundaryCondition = :open,STotZero = false)
    INDEX = 1:L
    if S==0.5
        JuMP.@variable(model, Sz[INDEX, INDEX], Bin)
    elseif S==1
        JuMP.@variable(model, -S<= Sz[INDEX, INDEX]<=S, Int)
    else
        error("S > 1 not implemented")
    end
    setConstraints!(model, L, opt, boundaryCondition,STotZero)
end

function setConstraints!(model, L, opt, boundaryCondition,STotZero)
    Sz = model[:Sz]
    for i = 1:L
        for j = 1:L
            iseven(i + j) || continue
            setConstraint!(model, Sz, i, j,opt,boundaryCondition,L)
        end
    end

    if STotZero
        JuMP.@constraint(model, sum(Sz) == 0)
    end
    return model
end

function setConstraint!(model, Sz, i, j,opt,boundaryCondition,L)

    if !plaquetteIsInBounds(Sz, i, j) &&boundaryCondition != :periodic
        return
    end

    neighborsites = (
        (i, j+1),
        (i-1, j+1),
        (i-1, j),
        (i-1, j-1),
        (i, j-1),
        (i+1, j-1),
        (i+1, j),
        (i+1, j+1),
    )
        
    s = map(neighborsites) do (i,j)
        CartesianIndex(wrap_indices_periodic(i,j,L,L,0))
    end
    # @info "" s
    # error("")
    JuMP.@constraint(
        model,
        Sz[s[1]] + Sz[s[2]] - Sz[s[3]] - Sz[s[4]] + Sz[s[5]] + Sz[s[6]] - Sz[s[7]] - Sz[s[8]] == 0
        # Sz[i, j+1] + Sz[i-1, j+1] - Sz[i-1, j] - Sz[i-1, j-1] + Sz[i, j-1] + Sz[i+1, j-1] -
        # Sz[i+1, j] - Sz[i+1, j+1] == 0
    )
    return nothing
end


function fixSpins!(model, fixInds, vals)

    Sz = model[:Sz]

    for i in eachindex(Sz)
        if JuMP.is_fixed(Sz[i])
            JuMP.unfix(Sz[i])
        end
    end
    for (I, v) in zip(fixInds, vals)
        JuMP.fix(Sz[I], v; force = true)
    end
    return model
end

function getRandomSpins(L::Integer, FixedFraction::Real;S=0.5)
    # numFixed = L^2 ÷ FixedFraction
    numFixed = round(Int, L^2 * FixedFraction)
    fixInds = unique!(rand(CartesianIndices((L, L)), numFixed))
    if S == 0.5
        vals = rand(0:1, length(fixInds))
    elseif S == 1
        vals = rand(-1:1, length(fixInds))
    else
        error("S > 1 not implemented")
    end
    return fixInds, vals
end

function setUpSpiderWeb(L::Integer,ENV;boundaryCondition=:open,S=0.5,STotZero=false, kwargs...)
    model = JuMP.Model(() -> Gurobi.Optimizer(ENV))
    JuMP.set_silent(model)
    setUpSpiderWeb!(model, L;boundaryCondition,S,STotZero)
    # JuMP.set_optimizer_attribute(model, "LogFile", LogFile)
    setGurobiParameters!(model; OutputFlag = false, TimeLimit = 3 * 60, kwargs...)
    return model
end


function setGurobiParameters!(model; kwargs...)
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
    for (k, v) in kwargs
        JuMP.set_optimizer_attribute(model, string(k), v)
    end

end

function constructGroundstate(model, L, DatType)
    Sz = model[:Sz]

    JuMP.optimize!(model)

    step_successful = JuMP.is_solved_and_feasible(model)
    term_status = JuMP.termination_status(model)
    if !step_successful
        return _default_arr_type(DatType,0,0), step_successful, term_status
    end

    sol = _default_arr_type(DatType,L,L)
    for i in eachindex(Sz)
        sol[i] = round(DatType, JuMP.value(Sz[i]))
    end

    return sol, step_successful, term_status
end

_default_arr_type(::Type{Bool},Lx,Ly) = BitMatrix(undef, Lx,Ly)
_default_arr_type(::Type{I},Lx,Ly) where I = Matrix{I}(undef, Lx,Ly)

function constructGroundstates(
    L::Integer,
    modelfactory::Function,
    NRuns::Integer,
    fixedFraction::Real,
    Dat_type::Type,
    SpinValue,
)
    p = ProgressMeter.Progress(NRuns; showspeed = true)

    generate_showvalues(iter, numSolutions) =
        () -> [(:iter, iter.value), (:numSolutions, numSolutions.value)]
    
    exampleSol = _default_arr_type(Dat_type,L,L)

    sols = Vector{typeof(exampleSol)}(undef, NRuns)
    # statuses = falses(NRuns)
    statuses = Vector{MOI.TerminationStatusCode}(undef,NRuns)
    numSolutions = Threads.Atomic{Int}(0)
    numfinished = Threads.Atomic{Int}(0)

    chunks = ChunkSplitters.chunks(1:NRuns,n=Threads.nthreads())
    Threads.@threads for inds in chunks
        ENV = Gurobi.Env()
        for iter in inds
            model = modelfactory(L,ENV)
            JuMP.set_attribute(model, MOI.NumberOfThreads(), 1)
            fixinds, vals = getRandomSpins(L, fixedFraction;S=SpinValue)
            fixSpins!(model, fixinds, vals)

            x, step_successful,status = constructGroundstate(model, L, Dat_type)
            sols[iter] = x
            statuses[iter] = status
            Threads.atomic_add!(numfinished, 1)
            if step_successful
                Threads.atomic_add!(numSolutions, 1)
            end
            ProgressMeter.next!(p; showvalues = generate_showvalues(numfinished, numSolutions))
        end
    end
    # filter!(!isempty, sols)
    filter!(!isempty, sols)
    return (;solutions = sols,statuses)

end

function constructGroundstates_noProgbar(
    L::Integer,
    modelfactory::Function,
    NRuns::Integer,
    fixedFraction::Real,
    Dat_type::Type,
    SpinValue,
)

    exampleSol = _default_arr_type(Dat_type,L,L)

    sols = Vector{typeof(exampleSol)}(undef, NRuns)
    # statuses = falses(NRuns)
    statuses = Vector{MOI.TerminationStatusCode}(undef,NRuns)
    numSolutions = Threads.Atomic{Int}(0)
    numfinished = Threads.Atomic{Int}(0)

    chunks = ChunkSplitters.chunks(1:NRuns,n=Threads.nthreads())
    Threads.@threads for inds in chunks
        ENV = Gurobi.Env()
        for iter in inds
            model = modelfactory(L,ENV)
            JuMP.set_attribute(model, MOI.NumberOfThreads(), 1)
            fixinds, vals = getRandomSpins(L, fixedFraction,S=SpinValue)
            fixSpins!(model, fixinds, vals)

            x, step_successful,status = constructGroundstate(model, L, Dat_type)
            sols[iter] = x
            statuses[iter] = status
            Threads.atomic_add!(numfinished, 1)
            if step_successful
                Threads.atomic_add!(numSolutions, 1)
            end
            ProgressMeter.next!(p; showvalues = generate_showvalues(numfinished, numSolutions))
        end
    end
    # filter!(!isempty, sols)
    filter!(!isempty, sols), 
    return (;solutions = sols,statuses)

end


function constructGroundstates(
    L::Integer,
    NRuns::Integer,
    fixedFraction::Real;
    TimeLimit = 20,
    boundaryCondition = :open,
    kwargs...,

)
    modelfactory(L,Env) = setUpSpiderWeb(L,Env; TimeLimit, OutputFlag = false, S=0.5,boundaryCondition,kwargs...)

    constructGroundstates(L, modelfactory, NRuns, fixedFraction::Real,Bool,0.5)
end

function constructGroundstatesSpin1(
    L::Integer,
    NRuns::Integer,
    fixedFraction::Real;
    TimeLimit = 20,
    boundaryCondition = :open,
    progress = true,
    kwargs...,

)
    # modelfactory(L,Env = GRB_ENV_REF[]) = setUpSpiderWeb(L; TimeLimit, OutputFlag = false, S=1,boundaryCondition,kwargs...)
    modelfactory(L,Env) = setUpSpiderWeb(L,Env; TimeLimit, OutputFlag = false, S=1,boundaryCondition,kwargs...)
    if progress
        return constructGroundstates(L, modelfactory, NRuns, fixedFraction,Int8,1)
    else
        return constructGroundstates_noProgbar(L, modelfactory, NRuns, fixedFraction,Int8,1)
    end
end