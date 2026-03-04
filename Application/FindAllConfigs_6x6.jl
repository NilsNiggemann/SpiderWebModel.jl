using SpiderWebModel
using SpiderWebModel.JuMP
using SpiderWebModel.Gurobi
import SpiderWebModel as SW
# Function to find all configurations for a 6x6 lattice using Gurobi
function find_all_configs(L)
    println("Setting up Gurobi problem to find all configurations for a $(L)x$(L) lattice...")

    # Create a Gurobi model
    model = JuMP.Model(Gurobi.Optimizer)
    
    # Configure Gurobi to find multiple solutions
    JuMP.set_optimizer_attribute(model, "PoolSearchMode", 2)  # Search for multiple solutions
    # JuMP.set_optimizer_attribute(model, "PoolSolutions", 2000000000)  # Store up to 2 billion solutions
    JuMP.set_optimizer_attribute(model, "PoolSolutions", 10000000)  # More reasonable
    
    # 1. Add a time limit (e.g., 1 hour)
    JuMP.set_optimizer_attribute(model, "TimeLimit", 3600)

    # 3. Add solution limit to stop early
    JuMP.set_optimizer_attribute(model, "SolutionLimit", 100_000_000)

    # 4. Use MIPFocus to find solutions faster (less emphasis on proving optimality)
    JuMP.set_optimizer_attribute(model, "MIPFocus", 1)

    SW.setUpSpiderWeb!(model, L; boundaryCondition=:periodic, S=1)
    JuMP.optimize!(model)

    # Check the status of the optimization
    if JuMP.termination_status(model) ∉ [MOI.OPTIMAL, MOI.SOLUTION_LIMIT]
        println("Optimization did not find solutions. Status: $(JuMP.termination_status(model))")
        return []
    end

    # Extract all solutions from the solution pool
    n_solutions = JuMP.result_count(model)
    println("Found $n_solutions solutions")
    
    solutions = []
    for sol_idx in 1:n_solutions
        Sz = model[:Sz]
        solution = [JuMP.value(Sz[i, j]; result=sol_idx) for i in 1:L, j in 1:L]
        push!(solutions, solution)
    end

    return solutions
end

# Run the function
solutions = find_all_configs(6)
println("Total configurations found: $(length(solutions))")