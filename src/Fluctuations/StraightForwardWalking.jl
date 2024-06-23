abstract type AbstractOperator end
operatorname(X::T) where T <: AbstractOperator = string(T)
struct PlaquetteFlipOperator <: AbstractOperator 
    AffectedPlaquettes::Matrix{OrderedCollections.OrderedSet{Int}}
end
operatorname(X::PlaquetteFlipOperator) = "PlaquetteFlipOperator"

function PlaquetteFlipOperator(S::StencilSpinConfig)
    AffectedPlaquettes = precomputeAffectedPlaquettes(S)
    return PlaquetteFlipOperator(AffectedPlaquettes)
end

function apply_operator!(Walker::SpiderWebWalker,O::PlaquetteFlipOperator,weightfunc::T,I) where T
    
    weights = getWeightList1Move!(Walker,O.AffectedPlaquettes,weightfunc,I)

    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = (+1,-1)[moveidx]
    w = sum(weights)
    
    applyPlaquette!(Walker.Config, I[1], I[2], move[1])
    return w
end

function getWeightList1Move!(Walker::SpiderWebWalker,AffectedPlaquettes,weightfunc::T,I) where T
    P⁺,P⁻ = P_applicable(Walker.Config,I)
    P⁺ || P⁻ || return SA[0.,0.]

    getNPlaq!(Walker)

    w⁺ = P⁺ ? getWeight1Move!(Walker,AffectedPlaquettes,weightfunc,I,+1) : 0.
    w⁻ = P⁻ ? getWeight1Move!(Walker,AffectedPlaquettes,weightfunc,I,-1) : 0.

    return SA[w⁺,w⁻]
end

function getWeight1Move!(Walker::SpiderWebWalker,AffectedPlaquettes,weightfunc::T,I,move) where T
    (;Config,n_x,n_x´) = Walker

    i1,i2 = I
    applyPlaquette!(Config, i1,i2, move)
    
    indices = AffectedPlaquettes[i1,i2]
    getNPlaq!(Walker,indices)

    N□ = getNPlaq_difference(n_x,n_x´,indices) 

    weight = weightfunc(N□)
    applyPlaquette!(Config, i1,i2, -move)

    return weight
end

"""Operator which flips two plaquettes (I,J), where I is assumed to be the reference plaquette (e.g.) at the origin. Used to sample the
<PᵢPⱼ + PᵢPⱼ† + Pᵢ†Pⱼ + Pᵢ†Pⱼ†> correlator.
"""
struct BBOperator <: AbstractOperator
    I::Tuple{Int,Int}
    AffectedPlaquettes::Matrix{Vector{Int}}
end

# getBBOperator(I::Tuple,AffectedPlaquettes::Vector) =  BBOperator!(copy(AffectedPlaquettes),I)
function BBOperator!(AffectedPlaquettes::AbstractMatrix,Itup::Tuple)
    I = CartesianIndex(Itup)
    affI = AffectedPlaquettes[I]
    
    for J in CartesianIndices(AffectedPlaquettes)
        if J !== I && isassigned(AffectedPlaquettes,J)
            unique!(append!(AffectedPlaquettes[J],affI))
        end
    end
    return BBOperator(Itup,AffectedPlaquettes)
end
function BBOperator(S::StencilSpinConfig,I::Tuple)
    AffectedPlaquettes = precomputeAffectedPlaquettes(S)
    return BBOperator!(AffectedPlaquettes,I)
end

function apply_operator!(Walker::SpiderWebWalker,O::BBOperator,weightfunc::T,J) where T
    I = O.I
    B2_flip_moves = SA[
        (1,1), #(+,+)
        (1,-1), #(+,-)
        (-1,1), #(-,+)
        (-1,-1), #(-,-)
    ]
    
    weights = getWeightList2Moves!(Walker,B2_flip_moves,O.AffectedPlaquettes,weightfunc,I,J)
    w = sum(weights)
    if iszero(w)
        return w
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))

    move = B2_flip_moves[moveidx]

    applyPlaquette!(Walker.Config, I[1], I[2], move[1])
    applyPlaquette!(Walker.Config, J[1], J[2], move[2])
    return w
end

function getWeightList2Moves!(Walker::SpiderWebWalker,moves,AffectedPlaquettes,weightfunc::T,I,J) where T
    getNPlaq!(Walker)

    weights = map(moves) do move
        getWeight2Moves!(Walker,AffectedPlaquettes,weightfunc,I,J,move)
    end
    return weights
end

function getWeight2Moves!(Walker::SpiderWebWalker,AffectedPlaquettes,weightfunc::T,I,J,move) where T
    (;Config,n_x,n_x´) = Walker

    idx_I,idx_J = @. 1 + (1-move) ÷ 2
    move_I, move_J = move

    P_applicable(Config,I)[idx_I] || return 0.
    i1,i2 = I
    applyPlaquette!(Config, i1,i2, move_I)
    
    if !P_applicable(Config,J)[idx_J] 
        applyPlaquette!(Config, i1,i2, -move_I)
        return 0.
    end
    j1,j2 = J
    applyPlaquette!(Config, j1,j2, move_J)
    
    indices = AffectedPlaquettes[j1,j2]
    getNPlaq!(Walker,indices)

    N□ = getNPlaq_difference(n_x,n_x´,indices) 

    weight = weightfunc(N□)
    applyPlaquette!(Config, j1,j2, -move_J)
    applyPlaquette!(Config, i1,i2, -move_I)

    return weight
end

function initialize_forward_walking!(Walkers,weights,O::AbstractOperator,Configs,J,weightfunc::T) where T
    # @inbounds for (α, Walker) in enumerate(Walkers)
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        ConfView = @view Configs[:,:,α]
        get_config(Walker) .= ConfView
        wa = apply_operator!(Walker,O,weightfunc,J)
        weights[α] = wa
    end
end

function straight_forward_walking!(setup,NSteps,nBranch,weightfunc::T,Λ,w_avg_estimate) where T
    
    (;AffectedPlaquetteList,Walkers,weights,TotalWeights,reconfiguration_buffer,reconfigurationTable) = setup

    Operator_weight = mean(weights)
    reconfigurationList = @view reconfigurationTable[:,1]
    reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    
    if all(iszero,weights)
        TotalWeights .= 0
        return TotalWeights
    end
    for i in 1:NSteps
        propagateWalkers!(Walkers,weights,AffectedPlaquetteList,weightfunc,Λ,nBranch,w_avg_estimate)

        TotalWeights[i] = mean(weights)
        reconfigurationList = @view reconfigurationTable[:,i]
        reconfiguration!(Walkers,reconfigurationList,reconfiguration_buffer,weights)
    end
    TotalWeights[begin] *= Operator_weight
    return TotalWeights
end

function setup_operatorObservables(mProj,NumObs,NSteps,Op::AbstractOperator,outfile::AbstractString)
    dataset_name = string("Weights",operatorname(Op))
    OperatorWeights = h5open(outfile,"cw") do file
        createMMapArray(file,dataset_name,Float64,(mProj,NumObs,NSteps))
    end
    return OperatorWeights
end

setup_operatorObservables(mProj,NumObs,NSteps,Op::AbstractOperator,outfile::Nothing) = zeros(NSteps,mProj,NumObs)

function measure_operator(InitialState,outfile,SaveConfigs,mProj,nBranch,O::AbstractOperator,weightfunc::T,Λ,w_avg_estimate,AllPlaqs = collect(plaquetteIterator(InitialState))) where T
    Lx,Ly,Nwalkers,NSteps = size(SaveConfigs)
    NSteps = NSteps
    setup = setup_many_walker_GFMC(InitialState,Nwalkers,mProj)
    results = setup_operatorObservables(mProj,length(AllPlaqs),NSteps,O,outfile)
    for n in 1:NSteps
        Configs = @view SaveConfigs[:,:,:,n]
        for (j,J) in enumerate(AllPlaqs)
            initialize_forward_walking!(setup.Walkers,setup.weights,O,Configs,J,weightfunc)
            res = straight_forward_walking!(setup,mProj,nBranch,weightfunc,Λ,w_avg_estimate)
            results[:,j,n] .= res

        end
    end
    return results
end

function measure_operator(InitialState,SaveConfigs,mProj,nBranch,O::AbstractOperator,weightfunc,Λ,w_avg_estimate,AllPlaqs = collect(plaquetteIterator(InitialState));outfile = nothing)
    measure_operator(InitialState,outfile,SaveConfigs,mProj,nBranch,O,weightfunc,Λ,w_avg_estimate,AllPlaqs)
end

function measure_operator(InitialState,inputfile::AbstractString,O::AbstractOperator,mProj::Integer,weightfunc,w_avg_estimate,AllPlaqs = collect(plaquetteIterator(InitialState));outfile = nothing)
    SaveConfigs = readMMapArray(inputfile,"SaveConfigs")
    nBranch = h5read(inputfile,"nBra")
    Λ = h5read(inputfile,"Λ")
    measure_operator(InitialState,SaveConfigs,mProj,nBranch,O,weightfunc,Λ,w_avg_estimate,AllPlaqs;outfile)
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