import SpiderWebModel as SW
using CairoMakie
##
Lx = 6
Ly = 7
@time a = SW.constructAllConfigs(Lx, Ly, SW.ALLGS_S12)
##
a_rec = SW.fillEmptyStates(a, Lx, Ly, SW.ALLGS_S12)
##
SW.plotApplPlaquettes(SW.reconstructTiling_xDirec(Lx, Ly, rand(a), SW.ALLGS_S12))
##
# SW.plotApplPlaquettes(SW.fillEmptyState(a_rec[3])[3])
##
Offset = 1
function getPeriodicTilings(states, offset)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(Lx, Ly, states[i], SW.ALLGS_S12, offset)
    end
    return findall(isValidTiling)
end

function getPeriodicTilings(states::AbstractVector{<:SW.SpinConfig}, offset)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(states[i], 5Lx, 5Ly, offset)
    end
    return findall(isValidTiling)
end

@time periodicTilings = getPeriodicTilings(a_rec, Offset)
##
function getDensity(state_lazy)
    state = SW.getPeriodicState_history(4Lx, 4Ly, Lx, Ly, state_lazy, SW.ALLGS_S12, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

function getDensity(state::SW.SpinConfig)
    state = SW.getPeriodicState(state, 70Lx, 70Ly, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end
densities = fetch.([Threads.@spawn getDensity(a_rec[x]) for x in periodicTilings])
count(!=(0), densities)

##
ind = periodicTilings[argmax(densities)]
Maxstate = SW.getPeriodicState(a_rec[ind], 4Lx, 4Ly, Offset)
SW.plotApplPlaquettes(Maxstate) |> display
##

##
maxdens = sort(densities, rev = true)[20]
inds = periodicTilings[findall(>=(maxdens), densities)]
let
    for ind in inds
        # densities[ind] == 0 && continue
        # state = SW.getPeriodicState_history(4Lx,4Ly,Lx,Ly,a[i],SW.ALLGS_S12)
        state = SW.getPeriodicState(a_rec[ind], 4Lx, 4Ly, Offset)
        SW.plotApplPlaquettes(state) |> display
    end
end
##
rddensities = round.(densities ./ 80, digits = 1)
uniquedens = sort!(unique(rddensities))
##
inds = [periodicTilings[findfirst(==(u), rddensities)] for u in uniquedens]

uniqueConfs = [SW.getPeriodicState(a_rec[ind], 3Lx, 3Ly, Offset) for ind in inds]

let
    N = length(inds)
    size = Int.((ceil(sqrt(N)), ceil(sqrt(N))))
    fig = Figure(size = 900 .* (size) ./ sqrt(N))

    exampleState = SW.getPeriodicState(a_rec[1], 3Lx, 3Ly, Offset)

    axes = [Axis(fig[i, j];
        SW.getConfigAxis(exampleState)...,
        xticklabelsvisible = i == size[2],
        yticklabelsvisible = j == 1) for i in 1:size[2], j in 1:size[1]]
    for (i, state) in enumerate(uniqueConfs)
        @assert SW.fulFillsConstraint(state)
        SW.plotApplPlaquettes!(axes[i], state)
        # SW.plotSpinConfig!(axes[i],state)
    end

    fig
end
