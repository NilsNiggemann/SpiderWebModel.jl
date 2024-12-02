struct PlaquetteFlipOperator{T} <: AbstractOperator 
    AffectedPlaquettes::Matrix{T}
end
operatorname(X::PlaquetteFlipOperator) = "PlaquetteFlipOperator"

function PlaquetteFlipOperator(S::StencilSpinConfig)
    AffectedPlaquettes = precomputeAffectedPlaquettes(S)
    return PlaquetteFlipOperator(AffectedPlaquettes)
end

function apply_operator!(Walker::SpiderWebWalker,O::PlaquetteFlipOperator,GuidingFuncBuffer,ψG::T,I) where T
    
    weights = getWeightList1Move!(Walker,GuidingFuncBuffer,ψG,I)

    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = (+1,-1)[moveidx]
    w = sum(weights)
    
    applyPlaquette!(Walker.Config, I[1], I[2], move[1])
    post_move_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker,move)
    return w
end

function getWeightList1Move!(Walker::SpiderWebWalker,GuidingFuncBuffer,ψG::T,I) where T
    i1,i2 = I
    P⁺,P⁻ = P_applicable(Walker.Config,i1,i2)
    P⁺ || P⁻ || return SA[0.,0.]

    getNPlaq!(Walker)

    w⁺ = P⁺ ? getWeight1Move!(Walker,GuidingFuncBuffer,ψG,(i1,i2,+1)) : 0.
    w⁻ = P⁻ ? getWeight1Move!(Walker,GuidingFuncBuffer,ψG,(i1,i2,-1)) : 0.

    return SA[w⁺,w⁻]
end

function getWeight1Move!(Walker::SpiderWebWalker,GuidingFuncBuffer,ψG::T,move) where T
    getNPlaq!(Walker)
    # weight = guidingfuncRatio(ψG,Walker,move,GuidingFuncBuffer)
    weight = guidingfuncRatio_naive(ψG,Walker,move)
    return weight
end

"""Operator which flips two plaquettes (I,J), where I is assumed to be the reference plaquette (e.g.) at the origin. Used to sample the
<PᵢPⱼ + PᵢPⱼ† + Pᵢ†Pⱼ + Pᵢ†Pⱼ†> correlator.
"""
struct BBOperator{T} <: AbstractOperator
    I::Tuple{Int,Int}
    AffectedPlaquettes::Matrix{T}
end

# getBBOperator(I::Tuple,AffectedPlaquettes::Vector) =  BBOperator!(copy(AffectedPlaquettes),I)
function BBOperator!(AffectedPlaquettes::AbstractMatrix,Itup::Tuple)
    I = CartesianIndex(Itup)
    affI = AffectedPlaquettes[I]
    
    for J in CartesianIndices(AffectedPlaquettes)
        if J !== I && isassigned(AffectedPlaquettes,J)
            union!(AffectedPlaquettes[J],affI)
        end
    end
    return BBOperator(Itup,AffectedPlaquettes)
end
function BBOperator(S::StencilSpinConfig,I::Tuple)
    AffectedPlaquettes = precomputeAffectedPlaquettes(S)
    return BBOperator!(AffectedPlaquettes,I)
end

function apply_operator!(Walker::SpiderWebWalker,O::BBOperator,GuidingFuncBuffer,ψG::T,J) where T
    I = O.I
    B2_flip_moves = SA[
        (1,1), #(+,+)
        (1,-1), #(+,-)
        (-1,1), #(-,+)
        (-1,-1), #(-,-)
    ]
    
    weights = getWeightList2Moves!(Walker,B2_flip_moves,GuidingFuncBuffer,ψG,I,J)
    w = sum(weights)
    if iszero(w)
        return w
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))

    move = B2_flip_moves[moveidx]

    move1 = (I[1],I[2],move[1])
    move2 = (J[1],J[2],move[2])
    applyPlaquette!(Walker.Config, move1)
    post_move_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker,move1)
    applyPlaquette!(Walker.Config, move2)
    post_move_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker,move2)

    return w
end

function getWeightList2Moves!(Walker::SpiderWebWalker,moves,GuidingFuncBuffer,ψG::T,I,J) where T
    premove_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker)

    weights = map(moves) do move
        getWeight2Moves!(Walker,GuidingFuncBuffer,ψG,I,J,move)
    end
    return weights
end

function getWeight2Moves!(Walker::SpiderWebWalker,GuidingFuncBuffer,ψG::T,I,J,move) where T
    (;Config) = Walker

    idx_I,idx_J = @. 1 + (1-move) ÷ 2
    move_I, move_J = move

    P_applicable(Config,I)[idx_I] || return 0.
    i1,i2 = I
    premove_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker)
    ψx = ψG(Walker)

    applyPlaquette!(Config, i1,i2, move_I)

    if !P_applicable(Config,J)[idx_J] 
        applyPlaquette!(Config, i1,i2, -move_I)
        return 0.
    end
    post_move_update_GWF_buffer!(GuidingFuncBuffer,ψG,Walker,(i1,i2,move_I))

    j1,j2 = J
    applyPlaquette!(Config, j1,j2, move_J)
    
    # indices = AffectedPlaquettes[j1,j2]

    # weight = guidingfuncRatio(ψG,n_x,n_x´,indices)
    ψx´ = ψG(Walker)
    weight = ψx´/ψx

    applyPlaquette!(Config, j1,j2, -move_J)
    applyPlaquette!(Config, i1,i2, -move_I)

    return weight
end

function initialize_forward_walking!(Walkers,weights,O::AbstractOperator,Configs,J,ψG::T,Guiding_function_buffer,nThreads) where T
    # @inbounds for (α, Walker) in enumerate(Walkers)
    batches = ChunkSplitters.chunks(eachindex(Walkers), n = nThreads)
    
    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        for α in αinds
            GWFBuffer = Guiding_function_buffer[i_chunk]
            Walker = Walkers[α]
            ConfView = @view Configs[:,:,α]
            get_config(Walker) .= ConfView
            compute_GWF_buffer!(GWFBuffer,ψG,Walker)
            wa = apply_operator!(Walker,O,GWFBuffer,ψG,J)
            weights[α] = wa
        end
    end
end
initialize_forward_walking!(Problem::AbstractGFMCProblem,O::AbstractOperator,Configs,OperatorList,nThreads) = initialize_forward_walking!(Problem.Walkers,Problem.weights,O,Configs,OperatorList,Problem.ψG,Problem.Guiding_function_buffer,nThreads)

function straight_forward_walking!(prob::AbstractGFMCProblem,TotalWeights,reconfigurationList,nThreads::Integer)
    
    (;Walkers,weights,Guiding_function_buffer,reconfiguration_buffer,ψG,method) = prob

    NSteps = size(TotalWeights,1)
    Operator_weight = mean(weights)

    reconfiguration!(Walkers,Guiding_function_buffer,reconfigurationList,reconfiguration_buffer,weights)
    if all(iszero,weights)
        TotalWeights .= 0
        return TotalWeights
    end

    for i in 1:NSteps

        propagateWalkers!(Walkers,weights,Guiding_function_buffer,nThreads,ψG,method)
        
        TotalWeights[i] = mean(weights)

        reconfiguration!(Walkers,Guiding_function_buffer,reconfigurationList,reconfiguration_buffer,weights)
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

setup_operatorObservables(mProj,NumObs,NSteps,Op::AbstractOperator,outfile::Nothing) = zeros(mProj,NumObs,NSteps)

function measure_operator(InitialState,method::AbstractGFMCMethod,outfile,SaveConfigs,mProj,O::AbstractOperator,ψG::T,AllPlaqs = collect(plaquetteIterator(InitialState)),nThreads = 2*Threads.nthreads()) where T
    Lx,Ly,Nwalkers,NSteps = size(SaveConfigs)

    setup = setup_many_walker_GFMC(InitialState,Nwalkers,nThreads)
    
    results = setup_operatorObservables(mProj,length(AllPlaqs),NSteps,O,outfile)

    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState,Nwalkers)

    Problem = SpiderwebGFMCProblem(method,InitialState,ψG,setup.Walkers,setup.weights,Guiding_function_buffer,setup.reconfiguration_buffer,results)
    
    reconfigurationList = zeros(Int,length(Problem.Walkers))
    for n in 1:NSteps
        Configs = @view SaveConfigs[:,:,:,n]
        for (j,J) in enumerate(AllPlaqs)
            initialize_forward_walking!(Problem,O,Configs,J,nThreads)
            
            TotalWeights = @view results[:,j,n]
            fill_all_Buffers!(Problem,nThreads)
            straight_forward_walking!(Problem,TotalWeights,reconfigurationList,nThreads)
            # results[:,j,n] .= res

        end
    end
    return results
end

function measure_operator(InitialState,method::AbstractGFMCMethod,SaveConfigs,mProj,O::AbstractOperator,ψG,AllPlaqs = collect(plaquetteIterator(InitialState));outfile = nothing,nThreads=2*Threads.nthreads())
    measure_operator(InitialState,method,outfile,SaveConfigs,mProj,O,ψG,AllPlaqs,nThreads)
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

struct RandomPlaquetteFlipOperator{T} <: AbstractOperator 
    AffectedPlaquettes::Matrix{T}
end
operatorname(X::RandomPlaquetteFlipOperator) = "RandomPlaquetteFlipOperator"

function RandomPlaquetteFlipOperator(S::StencilSpinConfig)
    AffectedPlaquettes = precomputeAffectedPlaquettes(S)
    return RandomPlaquetteFlipOperator(AffectedPlaquettes)
end

function apply_operator!(Walker::SpiderWebWalker,O::RandomPlaquetteFlipOperator,Guiding_function_buffer,ψG::T,::Any) where T
    compute_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
    weights = updateWeightList!(Walker,Guiding_function_buffer,ψG)
    moves = Walker.moves
    moveidx = StatsBase.sample(StatsBase.Weights(weights))
    move = moves[moveidx]
    w = sum(weights)
    applyPlaquette!(Walker.Config, move[1], move[2], move[3])
    return w
end

function measureObservables(S::StencilSpinConfig,BOp::AbstractOperator,OperatorList,p::Integer,Nwalkers::Integer,NSteps::Integer,method::AbstractGFMCMethod,ψG::AbstractGuidingFunction;nThreads=Threads.nthreads(), kwargs...)
    result = startManyWalkerGFMC(S,method,Nwalkers,NSteps,ψG;kwargs...)
    return measureObservables(result,S,BOp,OperatorList,p,method,ψG,nThreads)
end

function measureObservables(result::GFMCObservables,S::StencilSpinConfig,BOp::AbstractOperator,OperatorList,p::Integer,method::AbstractGFMCMethod,ψG::AbstractGuidingFunction,nThreads)
    resSFW = measure_operator(S,method,result.SaveConfigs,p,BOp,ψG,OperatorList;nThreads)
    Gnp = precomputeNormalizedAccWeight(result.TotalWeights,1,max(p,2))

    ObVal = stack([get_observables_sfw(Gnp,resSFW[:,j,:]',mean(result.TotalWeights)) for j in eachindex(IndexLinear(),OperatorList)])
    return ObVal
end