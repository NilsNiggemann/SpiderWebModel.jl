using SpiderWebModel
import SpiderWebModel as SW
using SpiderWebModel.JuMP
using SpiderWebModel.HiGHS
using MakieHelpers
using CairoMakie
##
spiderweb = Model(HiGHS.Optimizer)
set_silent(spiderweb)
Highs_resetGlobalScheduler(1)
set_attribute(model, MOI.NumberOfThreads(), 1)
##
Lx = 35
Ly = 35
unregister(spiderweb, :x)
@variable(spiderweb, x[i = 1:Lx, j = 1:Ly], Bin);
##
for i = 1:Lx  # For each row
    for j = 1:Ly  # and each column
        # Sum across all the possible digits. One and only one of the digits
        # can be in this cell, so the sum must be equal to one.
        iseven(i + j) || continue
        SW.plaquetteIsInBounds(x, i, j) || continue
        @constraint(
            spiderweb,
            x[i, j+1] + x[i-1, j+1] - x[i-1, j] - x[i-1, j-1] + x[i, j-1] + x[i+1, j-1] -
            x[i+1, j] - x[i+1, j+1] == 0
        )
    end
end

##
function getSols(N, spiderweb, x, Lx, Ly)
    success = falses(N)
    sols = [BitMatrix(zeros(Bool, Lx, Ly)) for _ = 1:N]

    # Threads.@threads 
    for iter = 1:N
        fixInds = rand(CartesianIndices(x), Lx * Ly ÷ 6)
        for i in eachindex(x)
            if is_fixed(x[i])
                unfix(x[i])
            end
        end

        for I in fixInds
            fix(x[I], rand(0:1); force = true)
        end
        @time optimize!(spiderweb)
        step_successful = is_solved_and_feasible(spiderweb)
        println(
            "step ",
            iter,
            " out of ",
            N,
            " progress: ",
            iter / N * 100,
            "%",
            " success: ",
            step_successful,
        )
        step_successful || continue
        success[iter] = true
        currentSol = sols[iter]
        for i in eachindex(currentSol)
            currentSol[i] = round(Bool, value(x[i]))
        end
    end
    deleteat!(sols, (!).(success))
    return sols
end

##
sols = getSols(10, spiderweb, x, Lx, Ly)
##
Confs = SW.floatSpinConfig.(SW.SpinConfig.(sols, 1 / 2))
SW.plotFractons(Confs[end])
##
Sq = SW.getEqualWeightStructureFac(Confs)
with_theme(theme_PiTicks()) do
    heatmap(Sq.k, Sq.k, real(Sq.Sq_k))

end
