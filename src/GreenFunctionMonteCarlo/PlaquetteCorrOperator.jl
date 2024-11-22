"""Operator Δ + Σᵢⱼ cos(q rᵢⱼ) <PᵢPⱼ + PᵢPⱼ† + Pᵢ†Pⱼ + Pᵢ†Pⱼ†>.
The delta is to make sure the the operator is positive definite.
"""
struct BBqOperator <: AbstractOperator 
end

operatorname(X::BBqOperator) = "BBqOperator"

function getAllTwoMoves(Walker::SpiderWebWalker)
    nx = Walker.n_x

    moves = copy(Walker.moves)
    # secondmoves = [empty(moves) for _ in eachindex(moves)]
    # secondmoves = empty!([moves])

    Config = get_config(Walker)
    getMoves!(moves,Config)
    newmoves = copy(moves)
    Allmoves = empty!([(DIAGONAL_MOVE_ID,DIAGONAL_MOVE_ID)])
    sizehint!(Allmoves, length(moves)^2)
    
    for move in moves
        i,j,s = move
        applyPlaquette!(Config, i, j, s)

        getMoves!(newmoves,Config)
        for secondmove in newmoves
            push!(Allmoves,(move,secondmove))
        end
        applyPlaquette!(Config, i, j, -s)
    end

    # return (;moves,secondmoves)
    return Allmoves
end

function apply_operator!(Walker::SpiderWebWalker,O::BBqOperator,Guiding_function_buffer,ψG::AbstractGuidingFunction,q::SVector)
    AllMoves = getAllTwoMoves(Walker)
    if isempty(AllMoves)
        return 0.
    end
    # weights = getWeightListAll2Moves!(Walker,ψG,O,q)
    weights = getWeightListAll2Moves!(Walker,AllMoves,ψG,Guiding_function_buffer,O,q)

    w = sum(weights)
    if iszero(w)
        return w
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights))

    move1,move2 = AllMoves[moveidx]

    applyPlaquette!(Walker.Config, move1...)
    applyPlaquette!(Walker.Config, move2...)
    return w
end

function apply_operator_buffer!(Walker::SpiderWebWalker,O::BBqOperator,Guiding_function_buffer,ψG::AbstractGuidingFunction,q::SVector,operatorBuffer)
    (;AllMoves,psiRatios,weights,Rij_x,Rij_y) = operatorBuffer
    if isempty(AllMoves)
        return 0.
    end

    qx,qy = q
    wsum = 0.
    
    LoopVectorization.@turbo for i in eachindex(psiRatios,weights,Rij_x,Rij_y)
        qr_ij = qx* Rij_x[i] + qy*Rij_y[i]
        weights[i] = psiRatios[i]*cos(qr_ij*0.5)^2
        wsum += weights[i]
    end

    if iszero(wsum)
        return wsum
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights,wsum))

    move1,move2 = AllMoves[moveidx]

    applyPlaquette!(Walker.Config, move1...)
    applyPlaquette!(Walker.Config, move2...)
    return wsum
end

function getWeightListAll2Moves!(Walker::SpiderWebWalker,AllMoves,ψG::AbstractGuidingFunction,Guiding_function_buffer,O::BBqOperator,q)
    # nx = getNPlaq!(Walker)
    compute_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
    weights = [getWeight2Moves!(Walker,ψG,move1,move2,O,q,Guiding_function_buffer) for (move1,move2) in AllMoves]

    return weights
end

function getWeight2Moves!(Walker::SpiderWebWalker,ψG::AbstractGuidingFunction,move1,move2,O::BBqOperator,q,Guiding_function_buffer)
    (;Config) = Walker

    premove_update_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
    ψx´_ψx = guidingfuncRatio(ψG,Walker,move1,Guiding_function_buffer) # ψx´/ψx

    applyPlaquette!(Config, move1)
    post_move_update_GWF_buffer!(Guiding_function_buffer,ψG,Walker,move1)

    premove_update_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
    ψx´´_ψx´ = guidingfuncRatio(ψG,Walker,move2,Guiding_function_buffer) # ψx´´/ψx´
    
    weight = ψx´´_ψx´ * ψx´_ψx

    i_x, i_y, move_I = move1
    j_x, j_y, move_J = move2

    qr_ij = q ⋅ SA[i_x-j_x, i_y-j_y]

    OperatorWeight = 2*(cos(qr_ij*0.5)^2)#/length(Walker.Plaquette_positions)

    post_move_update_GWF_buffer!(Guiding_function_buffer,ψG,Walker,inverse_move(move1)) # undo the first move
    applyPlaquette!(Config,inverse_move(move1))

    return weight*OperatorWeight
end

getRij(move1,move2) = SA[move1[1]-move2[1],move1[2]-move2[2]]

const MOVE_TYPE = Tuple{Int,Int,Int}

struct BBQ_WaveFunctionBuffer_1
    AllMoves::Vector{Tuple{MOVE_TYPE,MOVE_TYPE}}
    psiRatios::Vector{Float64}
    weights::Vector{Float64}
    Rij_x::Vector{Int}
    Rij_y::Vector{Int}
end

function buffer_BBQ_WFWeights(Walker::SpiderWebWalker,ψG::AbstractGuidingFunction,Guiding_function_buffer)
    AllMoves = getAllTwoMoves(Walker)
    psiRatios = getWeightListAll2Moves!(Walker,AllMoves,ψG,Guiding_function_buffer,BBqOperator(),SA[0.,0.])
    weights = similar(psiRatios) .= 0. 
    Rij = (getRij(move1,move2) for (move1,move2) in AllMoves)
    Rij_x = [rij[1] for rij in Rij]
    Rij_y = [rij[2] for rij in Rij]
    return BBQ_WaveFunctionBuffer_1(AllMoves,psiRatios,weights,Rij_x,Rij_y)
end

function buffer_BBQ_WFWeights(Walkers::Vector{<:SpiderWebWalker},ψG::AbstractGuidingFunction,Guiding_function_buffer)
    chunks = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())
    allbuffers = Vector{BBQ_WaveFunctionBuffer_1}(undef,length(Walkers))
    Threads.@threads for (i_chunk,αinds) in enumerate(chunks)
        GWFBuffer = Guiding_function_buffer[i_chunk]
        for α in αinds
            Walker = Walkers[α]
            allbuffers[α] = buffer_BBQ_WFWeights(Walker,ψG,Guiding_function_buffer[i_chunk])
        end
    end
    return allbuffers    
end

function measure_operator(InitialState,method::AbstractGFMCMethod,outfile,SaveConfigs,mProj,O::BBqOperator,ψG::T,Allqs) where T
    Lx,Ly,Nwalkers,NSteps = size(SaveConfigs)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers)
    
    results = setup_operatorObservables(mProj,length(Allqs),NSteps,O,outfile)

    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState)

    Problem = SpiderwebGFMCProblem(method,InitialState,ψG,setup.Walkers,setup.weights,Guiding_function_buffer,setup.reconfiguration_buffer,results)

    reconfigurationList = zeros(Int,length(Problem.Walkers))
    for n in 1:NSteps
        Configs = @view SaveConfigs[:,:,:,n]
        fillWalkers!(Problem.Walkers,Configs)
        WF_buffers = buffer_BBQ_WFWeights(Problem.Walkers,ψG,Guiding_function_buffer)
        for (j,q) in enumerate(Allqs)
            initialize_buffered_forward_walking!(Problem,O,Configs,q,WF_buffers)

            TotalWeights = @view results[:,j,n]
            straight_forward_walking!(Problem,TotalWeights,reconfigurationList)
            # results[:,j,n] .= res

        end
    end
    return results
end
function fillWalkers!(Walkers,Configs)
    Threads.@threads for α in eachindex(Walkers)
        Walker = Walkers[α]
        ConfView = @view Configs[:,:,α]
        get_config(Walker) .= ConfView
    end
end

function _initialize_buffered_forward_walking!(Walkers,weights,O::BBqOperator,Configs,q,ψG::T,Guiding_function_buffer,WF_buffers) where T
    # @inbounds for (α, Walker) in enumerate(Walkers)
    batches = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())
    
    # for (i_chunk,αinds) in enumerate(batches)
    Threads.@threads for (i_chunk,αinds) in enumerate(batches)
        for α in αinds
            operatorBuffer = WF_buffers[α]
            # operatorBuffer = buffer_BBQ_WFWeights(Walkers[α],ψG)
            GWFBuffer = Guiding_function_buffer[i_chunk]
            Walker = Walkers[α]
            ConfView = @view Configs[:,:,α]
            get_config(Walker) .= ConfView
            wa = apply_operator_buffer!(Walker,O,GWFBuffer,ψG,q,operatorBuffer)
            weights[α] = wa
        end
    end
end
initialize_buffered_forward_walking!(Problem::AbstractGFMCProblem,O::BBqOperator,Configs,operator,WF_buffers) = _initialize_buffered_forward_walking!(Problem.Walkers,Problem.weights,O,Configs,operator,Problem.ψG,Problem.Guiding_function_buffer,WF_buffers)
