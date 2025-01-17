import SpiderWebModel as SW
using CairoMakie
function getPeriodicTilings(states, offset)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(Lx, Ly, states[i], SW.ALLGS_S12, offset)
    end
    return findall(isValidTiling)
end

function getPeriodicTilings(states::AbstractVector{<:SW.SpinConfig}, Lx,Ly,offset)
    isValidTiling = zeros(Bool, length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(states[i], 5Lx, 5Ly, offset)
    end
    return findall(isValidTiling)
end
function getDensity(state_lazy,Lx,Ly,Offset)
    state = SW.getPeriodicState_history(4Lx, 4Ly, Lx, Ly, state_lazy, SW.ALLGS_S12, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

function getDensity(state::SW.SpinConfig,Lx,Ly,Offset)
    state = SW.getPeriodicState(state, 70Lx, 70Ly, Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

##
Lx = 10
Ly = 10
@time a = SW.constructAllConfigs(Lx, Ly, SW.ALLGS_S12)
##
@time a_rec = SW.fillEmptyStates(a, Lx, Ly, SW.ALLGS_S12)
##
SW.plotApplPlaquettes(SW.reconstructTiling_xDirec(Lx, Ly, rand(a), SW.ALLGS_S12))
##
# SW.plotApplPlaquettes(SW.fillEmptyState(a_rec[3])[3])
##
Offset = 3


@time periodicTilings = getPeriodicTilings(a_rec, Offset)
##
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

    axes = [
        Axis(
            fig[i, j];
            SW.getConfigAxis(exampleState)...,
            xticklabelsvisible = i == size[2],
            yticklabelsvisible = j == 1,
        ) for i = 1:size[2], j = 1:size[1]
    ]
    for (i, state) in enumerate(uniqueConfs)
        @assert SW.fulFillsConstraint(state)
        SW.plotApplPlaquettes!(axes[i], state)
        # SW.plotSpinConfig!(axes[i],state)
    end

    fig
end
##
function ConstructInequivPeriodicConfigs(maxL)
    uniqueUCS = empty!(
        Dict(SW.PeriodicMatrix(zeros(Int8,maxL,maxL),2,2,0) => 1.)
    )

    for Lx in 1:maxL
        for Ly in 1:maxL
            @info "" Lx Ly length(uniqueUCS)
            @time "constructing Configs" UCpth = SW.constructAllConfigs(Lx, Ly, SW.ALLGS_S12)
            @time "fill empty sites" UC_Full = SW.fillEmptyStates(UCpth, Lx, Ly, SW.ALLGS_S12)

            # @time UCs = [Int8.(2 .*x) for x in SW.fillEmptyStates(a, Lx, Ly, SW.ALLGS_S12)]

            for Offset in 0:max(Lx,Ly)
                @time "Offset = $Offset" periodicTilings = getPeriodicTilings(UC_Full,Lx,Ly, Offset)
                @time "compute densities" densities = fetch.([Threads.@spawn getDensity(UC_Full[x],20,20,Offset) for x in periodicTilings])
                rddensities = round.(densities ./ 80, digits = 1)
                uniquedens = sort!(unique(rddensities))

                for u in uniquedens
                    u == 0 && continue
                    ind = periodicTilings[findfirst(==(u), rddensities)]

                    UnConf = SW.PeriodicMatrix(Int8.( 2 .* UC_Full[ind]), Lx,Ly, Offset)
                    uniqueUCS[UnConf] = u
                end
            end
        end
    end
    return uniqueUCS

end

uniqueUCS = ConstructInequivPeriodicConfigs(6)

##
# function findEnlargedUC(permat)
#     Lx,Ly = size(permat.UC)
    
#     newUC = permat[1:Lx*(permat.offset+1),1:Ly]
# end
# enlargedUCS = findEnlargedUC.(uniqueUCS)
# function getDens(UC)
#     S = SW.stencilConfig(0.5*UC, 1 / 2,boundaryCondition = :periodic)
#     return length(SW.getApplicablePlaquettes(S))
# end
# densities = getDens.(enlargedUCS)
# perm = sortperm(densities,rev=true)

# uniqueUCS = uniqueUCS[perm]
# densities = densities[perm]
# enlargedUCS = enlargedUCS[perm]


ds = collect(values(uniqueUCS))
perm = sortperm(ds,rev=true)
confs = collect(keys(uniqueUCS))[perm]

##
let
    N = length(uniqueUCS)
    size = Int.((ceil(sqrt(N)), ceil(sqrt(N))))
    fig = Figure(size = 900 .* (size) ./ sqrt(N))

    S = SW.stencilConfig(0.5*ones(12,12), 1 / 2)

    axs = [
        Axis(
            fig[i, j];
            SW.getConfigAxis(S)...,
            xticklabelsvisible = i == size[2],
            yticklabelsvisible = j == 1,
        ) for i = 1:size[2], j = 1:size[1]
    ]

    for (i, state) in enumerate(confs)
        S .= state[axes(S,1),axes(S,2)]
        # @assert SW.fulFillsConstraint(S)
        SW.plotApplPlaquettes!(axs[i], S)
        # SW.plotSpinConfig!(axes[i],state)
    end

    fig
end
