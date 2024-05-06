function get_B2_moves!(Walker::SpiderWebWalker,I)
    moves = Walker.moves
    empty!(moves)
    P1,P2 = P_applicable(Walker.Config,I)
    P1 && push!(moves,(I[1],I[2],1))
    P2 && push!(moves,(I[1],I[2],-1))
    return moves
end

function apply_B2_Operator!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc,I,J)
    w = 1.
    for site in (I,J)
        get_B2_moves!(Walker,site)
        if isempty(Walker.moves) 
            empty!(Walker.weights)
            push!(Walker.weights,0)
            return 0. #(0,0,0),Walker.weights
        end
        weights = getWeightList!(Walker,AffectedPlaquetteList,weightfunc,0)
        w *= sum(weights)
        moveidx = StatsBase.sample(StatsBase.Weights(weights))
        move = Walker.moves[moveidx]
        applyPlaquette!(Walker.Config, move[1], move[2], move[3])
    end
    return w
end

function initialize_forward_walkingB2!(setup,Configs,I,J,weightfunc)
    (;AffectedPlaquetteList,Walkers,weights) = setup

    for (α, Walker) in enumerate(Walkers)
        get_config(Walker) .= @view Configs[:,:,α]
        weights[α] = apply_B2_Operator!(Walker,AffectedPlaquetteList,weightfunc,I,J)
    end
    return setup
end

function straight_forward_walkingB2!(setup,NSteps,nBranch,Configs,weightfunc::T,I,J,Λ) where T
    
    initialize_forward_walkingB2!(setup,Configs,I,J,weightfunc)
    
    (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationList) = setup

    Operator_weight = mean(weights)
    reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    
    if all(iszero,weights)
        TotalWeights .= 0
        return TotalWeights
    end

    for i in 1:NSteps
        propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc,Λ,nBranch)

        TotalWeights[i] = mean(weights)
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    end
    TotalWeights[begin] *= Operator_weight
    return TotalWeights
end

function measure_B2_correlators(InitialState,saveConfigs,mProj,nBranch,I,weightfunc,Λ)
    Lx,Ly,Nwalkers,NSteps = size(saveConfigs)
    NSteps = NSteps
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,mProj)
    AllPlaqs = plaquetteIterator(first(setup.Walkers))
    results = zeros(NSteps,mProj,length(AllPlaqs))
    for n in 1:NSteps
        Configs = @view saveConfigs[:,:,:,n]
        for (j,J) in enumerate(AllPlaqs)
            # J == (5,6) || continue
            res = straight_forward_walkingB2!(setup,mProj,nBranch,Configs,weightfunc,I,J,Λ)
            results[n,:,j] .= res

        end
    end
    return results
end

function get_observables_sfw(Gnp,sfw_weights,meanweight,m = size(sfw_weights,2))
    num = zeros(m)
    denom = zeros(m)
    for p in 1:m
        for n in p:lastindex(Gnp,1)
            num[p] += prod(sfw_weights[n,i]/meanweight for i in 1:p)
            denom[p] += Gnp[n,p]
        end
    end
    return num./denom
end