function apply_B2_Operator!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::T,I,J) where T
    moves,weights = getWeightList2Moves!(Walker,AffectedPlaquetteList,weightfunc,I,J)
    w = sum(weights)
    if iszero(w)
        return w
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))

    move = moves[moveidx]

    applyPlaquette!(Walker.Config, I[1], I[2], move[1])
    applyPlaquette!(Walker.Config, J[1], J[2], move[2])
    return w
end

function modify_AffectedPlaquettes_B2!(AffectedPlaquettes,Itup::Tuple{Int,Int})
    I = CartesianIndex(Itup)
    affI = AffectedPlaquettes[I]
    
    for J in CartesianIndices(AffectedPlaquettes)
        if J !== I
            unique!(append!(AffectedPlaquettes[J],affI))
            # AffectedPlaquettes[J] = affI ∪ AffectedPlaquettes[J]
        end
    end
    return AffectedPlaquettes
end

modify_AffectedPlaquettes_B2(AffectedPlaquettes,I) =  modify_AffectedPlaquettes_B2!(copy(AffectedPlaquettes),I)

function getWeightList2Moves!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::T,I,J) where T
    getNPlaq!(Walker)

    B2_flip_moves = SA[
        (1,1), #(+,+)
        (1,-1), #(+,-)
        (-1,1), #(-,+)
        (-1,-1), #(-,-)
    ]
    weights = map(B2_flip_moves) do move
        getWeight2Moves!(Walker,AffectedPlaquetteList,weightfunc,I,J,move)
    end
    return B2_flip_moves, weights
end

function getWeight2Moves!(Walker::SpiderWebWalker,AffectedPlaquetteList,weightfunc::T,I,J,move) where T
    (;Config,n_x,n_x´) = Walker

    idx_I,idx_J = @. 1 + (1-move) ÷ 2
    move_I, move_J = move

    P_applicable(Config,I)[idx_I] || return 0.
    i1,i2 = I
    applyPlaquette!(Config, i1,i2, move[1])
    
    if !P_applicable(Config,J)[idx_J] 
        applyPlaquette!(Config, i1,i2, -move[1])
        return 0.
    end
    j1,j2 = J
    applyPlaquette!(Config, j1,j2, move[2])
    
    indices = AffectedPlaquetteList[j1,j2]
    getNPlaq!(Walker,indices)

    N□ = getNPlaq_difference(n_x,n_x´,indices) 

    weight = weightfunc(N□)
    applyPlaquette!(Config, j1,j2, -move[2])
    applyPlaquette!(Config, i1,i2, -move[1])

    return weight
end

function initialize_forward_walkingB2!(Walkers,weights,AffectedPlaquetteListB2,Configs,I,J,weightfunc::T) where T
    # @inbounds for (α, Walker) in enumerate(Walkers)
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        ConfView = @view Configs[:,:,α]
        get_config(Walker) .= ConfView
        wa = apply_B2_Operator!(Walker,AffectedPlaquetteListB2,weightfunc,I,J)
        weights[α] = wa
    end
end

function straight_forward_walkingB2!(setup,NSteps,nBranch,weightfunc::T,Λ) where T
    
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

function measure_B2_correlators(InitialState,saveConfigs,mProj,nBranch,I,weightfunc::T,Λ) where T
    Lx,Ly,Nwalkers,NSteps = size(saveConfigs)
    NSteps = NSteps
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,mProj)
    AffectedPlaquetteListB2 = modify_AffectedPlaquettes_B2(setup.AffectedPlaquetteList,I)

    AllPlaqs = plaquetteIterator(first(setup.Walkers))
    results = zeros(NSteps,mProj,length(AllPlaqs))
    for n in 1:NSteps
        Configs = @view saveConfigs[:,:,:,n]
        for (j,J) in enumerate(AllPlaqs)
            initialize_forward_walkingB2!(setup.Walkers,setup.weights,AffectedPlaquetteListB2,Configs,I,J,weightfunc)
            res = straight_forward_walkingB2!(setup,mProj,nBranch,weightfunc,Λ)
            results[n,:,j] .= res

        end
    end
    return results
end

function get_observables_sfw(Gnp,sfw_weights,meanweight,m = size(sfw_weights,2))
    num = zeros(m)
    denom = zeros(m)
    Gnp´ = cumprod(sfw_weights./meanweight,dims=2)
    for p in 1:m
        for n in p:lastindex(Gnp,1)
            num[p] += Gnp´[n,p]
            denom[p] += Gnp[n,p]
        end
    end
    return num./denom
end