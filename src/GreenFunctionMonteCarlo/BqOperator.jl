"""Operator Σᵢ cos(q rᵢ) <Pᵢ+ Pᵢ†>.
"""
struct BqOperator <: AbstractOperator 
end

operatorname(X::BqOperator) = "BqOperator"

struct BQ_WaveFunctionBuffer
    AllMoves::Vector{MOVE_TYPE}
    psiRatios::Vector{Float64}
    weights::Vector{Float64}
    Ri_x::Vector{Int}
    Ri_y::Vector{Int}
end

function apply_operator_buffer!(Walker::SpiderWebWalker,O::BqOperator,Guiding_function_buffer,ψG::AbstractGuidingFunction,q::SVector,operatorBuffer)
    (;AllMoves,psiRatios,weights,Ri_x,Ri_y) = operatorBuffer

    if isempty(AllMoves)
        return 0.
    end
    qx,qy = q
    wsum = 0.
    LoopVectorization.@turbo for i in eachindex(psiRatios,weights,Ri_x,Ri_y)
        qr_i = qx* Ri_x[i] + qy*Ri_y[i]
        weights[i] = 2*psiRatios[i]*cos(qr_i*0.5)^2
        wsum += weights[i]
    end
    if iszero(wsum)
        return wsum
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights,wsum))

    move1 = AllMoves[moveidx]

    applyPlaquette!(Walker.Config, move1...)
    return wsum
end

function buffer_BQ_WFWeights(Walker::SpiderWebWalker,ψG::AbstractGuidingFunction,Guiding_function_buffer)
    AllMoves = copy(getMoves!(Walker))
    compute_GWF_buffer!(Guiding_function_buffer,ψG,Walker)
    psiRatios = copy(updateWeightList!(Walker,Guiding_function_buffer,ψG))

    weights = similar(psiRatios) .= 0. 
    Ri_x = [ri[1] for ri in AllMoves]
    Ri_y = [ri[2] for ri in AllMoves]
    return BQ_WaveFunctionBuffer(AllMoves,psiRatios,weights,Ri_x,Ri_y)
end

function buffer_BQ_WFWeights(Walkers::Vector{<:SpiderWebWalker},ψG::AbstractGuidingFunction,Guiding_function_buffer)
    chunks = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())
    allbuffers = Vector{BQ_WaveFunctionBuffer}(undef,length(Walkers))
    Threads.@threads for (i_chunk,αinds) in enumerate(chunks)
        GWFBuffer = Guiding_function_buffer[i_chunk]
        for α in αinds
            Walker = Walkers[α]
            allbuffers[α] = buffer_BQ_WFWeights(Walker,ψG,GWFBuffer)
        end
    end
    return allbuffers    
end

function measure_operator(InitialState,method::AbstractGFMCMethod,outfile,SaveConfigs,mProj,O::BqOperator,ψG::T,Allqs) where T
    Lx,Ly,Nwalkers,NSteps = size(SaveConfigs)
    setup = setup_many_walker_GFMC(InitialState,Nwalkers)
    
    results = setup_operatorObservables(mProj,length(Allqs),NSteps,O,outfile)

    Guiding_function_buffer = allocate_GWF_buffers_threads(ψG,InitialState)
    Problem = SpiderwebGFMCProblem(method,InitialState,ψG,setup.Walkers,setup.weights,Guiding_function_buffer,setup.reconfiguration_buffer,results)

    reconfigurationList = zeros(Int,length(Problem.Walkers))
    for n in 1:NSteps
        Configs = @view SaveConfigs[:,:,:,n]
        fillWalkers!(Problem.Walkers,Configs)
        WF_buffers = buffer_BQ_WFWeights(Problem.Walkers,ψG,Guiding_function_buffer)
        for (j,q) in enumerate(Allqs)

            initialize_buffered_forward_walking!(Problem,O,Configs,q,WF_buffers)
            TotalWeights = @view results[:,j,n]
            straight_forward_walking!(Problem,TotalWeights,reconfigurationList)
            # results[:,j,n] .= res
        end
    end
    return results
end
