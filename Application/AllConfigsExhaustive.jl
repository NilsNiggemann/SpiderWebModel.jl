import SpiderWebModel as SW
using CairoMakie
##
Lx = 7
Ly = 6
@time a = SW.constructAllConfigs(Lx,Ly,SW.ALLGS_S12)
##
a_rec = SW.fillEmptyStates(a,Lx,Ly,SW.ALLGS_S12)
##
SW.plotApplPlaquettes(SW.reconstructTiling_xDirec(Lx,Ly,rand(a),SW.ALLGS_S12))
##
# SW.plotApplPlaquettes(SW.fillEmptyState(a_rec[3])[3])
##
Offset = 1
function getPeriodicTilings(states,offset)
    isValidTiling = zeros(Bool,length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(Lx,Ly,states[i],SW.ALLGS_S12,offset)
    end
    return findall(isValidTiling)
end

function getPeriodicTilings(states::AbstractVector{<:SW.SpinConfig},offset)
    isValidTiling = zeros(Bool,length(states))
    Threads.@threads for i in eachindex(states)
        isValidTiling[i] = SW.isPeriodicTiling(states[i],5Lx,5Ly,offset)
    end
    return findall(isValidTiling)
end

@time periodicTilings = getPeriodicTilings(a_rec,Offset)
##
function getDensity(state_lazy)
    state = SW.getPeriodicState_history(4Lx,4Ly,Lx,Ly,state_lazy,SW.ALLGS_S12,Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end

function getDensity(state::SW.SpinConfig)
    state = SW.getPeriodicState(state,4Lx,4Ly,Offset)
    density = SW.getApplicablePlaquettes(state) |> length
end
densities = [getDensity(a_rec[x]) for x in periodicTilings]
count(!=(0),densities)

##
ind = periodicTilings[argmax(densities)]
Maxstate = SW.getPeriodicState(a_rec[ind],4Lx,4Ly,Offset)
SW.plotApplPlaquettes(Maxstate) |> display
##

##
let 
    inds = sort(densities,rev=true)[1:10]
    for ind in inds
        densities[ind] == 0 && continue
        # state = SW.getPeriodicState_history(4Lx,4Ly,Lx,Ly,a[i],SW.ALLGS_S12)
        state = SW.getPeriodicState(a_rec[ind],4Lx,4Ly,Offset)
        SW.plotApplPlaquettes(state) |> display
    end
end