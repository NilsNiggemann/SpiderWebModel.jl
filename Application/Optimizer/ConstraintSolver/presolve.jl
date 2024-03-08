using JuMP, Gurobi
include("../SpiderwebDefs.jl")
##
jump_model = direct_model(Gurobi.Optimizer())
setUpSpiderWeb!(jump_model,40)
setup = getRandomSpins(L,6)
fixSpins!(jump_model,setup...)
##
JuMP.set_optimizer_attribute(jump_model, "Presolve", 2) #aggressive
JuMP.set_optimizer_attribute(jump_model, "PrePasses", 2_000_000_000) #many prepasses
JuMP.set_optimizer_attribute(jump_model, "Aggregate", 2) #aggressive
JuMP.set_optimizer_attribute(jump_model, "AggFill", 2_000_000_000) #aggressive
JuMP.set_optimizer_attribute(jump_model, "PreSparsify", 2) #aggressive
##
JuMP.set_optimizer_attribute(jump_model, "Presolve", 0) #aggressive
JuMP.set_optimizer_attribute(jump_model, "PrePasses", 0) #many prepasses
JuMP.set_optimizer_attribute(jump_model, "Aggregate", 0) #aggressive
JuMP.set_optimizer_attribute(jump_model, "AggFill", 0) #aggressive
JuMP.set_optimizer_attribute(jump_model, "PreSparsify", 0) #aggressive
##
presolvedP = Ref{Ptr{Cvoid}}(C_NULL)
GRBpresolvemodel(backend(jump_model), presolvedP)
GRBwrite(presolvedP[], "presolve.lp")

##
model = backend(jump_model)
# GRBtunemodel(model)
# Gurobi.build(model)  # this instantiates the Gurobi model prior to solve
# solve(model)