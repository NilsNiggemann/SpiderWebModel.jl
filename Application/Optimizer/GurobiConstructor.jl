using SpiderWebModel
import SpiderWebModel as SW
using SpiderWebModel.JuMP
using SpiderWebModel.HiGHS
using MakieHelpers
using CairoMakie
using SpiderWebModel.Gurobi
using ProgressMeter
##
function setUpSpiderWeb(
    L;
    OutputFlag = false,
    TimeLimit = 3 * 60,
    LogFile = "ConfigSearch_$L.log",
)
    model = Model(Gurobi.Optimizer)
    Lx = Ly = L
    INDEX = 1:L
    @variable(model, y[INDEX, INDEX] >= 0, Bin)

    for i = 1:Lx
        for j = 1:Ly
            iseven(i + j) || continue
            SW.plaquetteIsInBounds(y, i, j) || continue
            @constraint(
                model,
                y[i, j+1] + y[i-1, j+1] - y[i-1, j] - y[i-1, j-1] +
                y[i, j-1] +
                y[i+1, j-1] - y[i+1, j] - y[i+1, j+1] == 0
            )
        end
    end
    ## Gurobi parameters - see https://www.gurobi.com/documentation/9.0/refman/finding_multiple_solutions.html
    JuMP.set_optimizer_attribute(model, "PoolSearchMode", 0)   # 1: non-exhaustive search mode
    JuMP.set_optimizer_attribute(model, "OutputFlag", OutputFlag)  # 0: suppress output
    JuMP.set_optimizer_attribute(model, "TimeLimit", TimeLimit)  # time limit before terminating in seconds
    if OutputFlag
        JuMP.set_optimizer_attribute(model, "LogFile", LogFile)
        JuMP.set_optimizer_attribute(model, "LogToConsole", 0)  # number of solutions

    end

    return model
end

function demoInitializer(L, initializeDenominator = 12)
    Conf = SW.SpinConfig(fill(NaN, L, L), 1 / 2)
    numFixed = L^2 ÷ initializeDenominator
    fixInds = rand(CartesianIndices(Conf), numFixed)
    vals = rand(-0.5:0.5, numFixed)
    for (I, v) in zip(fixInds, vals)
        Conf[I] = v
    end
    SW.plotSpinConfig(Conf)
end

function constructSW(model, L, numSols, initializeDenominator::Int)
    y = model[:y]
    numFixed = L^2 ÷ initializeDenominator
    fixInds = rand(CartesianIndices(y), numFixed)
    vals = rand(0:1, numFixed)
    INDEX = 1:L
    for i in eachindex(y)
        if is_fixed(y[i])
            unfix(y[i])
        end
    end
    for (I, v) in zip(fixInds, vals)
        fix(y[I], v; force = true)
    end

    JuMP.set_optimizer_attribute(model, "PoolSolutions", numSols)
    JuMP.optimize!(model)

    num_results = result_count(model)

    Results = BitMatrix[]
    step_successful = is_solved_and_feasible(model)
    if !step_successful
        return Results
    end
    for n = 1:num_results
        sol = BitMatrix(undef, L, L)
        for i in INDEX
            for j in INDEX
                sol[i, j] = round(Bool, JuMP.value(y[i, j]; result = n))
            end
        end
        push!(Results, sol)
    end

    return Results
end

function multipleRuns(model, L, NTot, NIndiv, initializeDenominator = 12)
    sols = BitMatrix[]
    p = Progress(NTot)
    # x = [BitMatrix(undef,0,0)]
    generate_showvalues(iter, x) = () -> [(:iter, iter), (:numSolutions, length(x))]

    @showprogress dt = 0.1 desc = "obtain configs..." for iter = 1:NTot
        x = constructSW(model, L, NIndiv, initializeDenominator)
        # @info "" i length(x)
        sols = vcat(sols, x)
        next!(p; showvalues = generate_showvalues(iter, sols))
    end
    return sols
end
##
L = 40
model =
    setUpSpiderWeb(L, OutputFlag = true, TimeLimit = 15, LogFile = "ConfigSearch_$L.log")
sols = multipleRuns(model, L, 50, 1, 7)
##
oldsols = h5read("Configs_$L.h5", "Confs")
newsols = cat(oldsols, Array{Bool,3}(stack(sols, dims = 3)), dims = 3)
##
# rm("Configs_$L.h5")
# h5write("Configs_$L.h5", "Confs", newsols)
##
Confs = SW.floatSpinConfig.(SW.SpinConfig.(eachslice(newsols, dims = 3), 1 / 2))
##
SW.plotFractons(rand(Confs))
##
Sq = SW.getEqualWeightStructureFac(Confs)
with_theme(theme_PiTicks()) do
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)

    hm = heatmap!(ax, Sq.k, Sq.k, real(Sq.Sq_k))
    @info "" sumrule = sum(real(Sq.Sq_k[1:end-1, 1:end-1])) sumRuleExpected = L^2 / 4
    Colorbar(fig[1, 2], hm)
    fig
    save("Application/figs/Sq_$L.pdf", fig)
end
