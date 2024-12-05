"""Operator Σᵢ cos(q rᵢ) <Pᵢ+ Pᵢ†>.
"""
struct BqOperator <: AbstractOperator 
    phase::Float64
end
BqOperator() = BqOperator(0.)

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
        weights[i] = 2*psiRatios[i]*cos(qr_i*0.5+O.phase)^2
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

    I0 = getCentralPlaquette(get_config(Walker))
    weights = similar(psiRatios) .= 0. 
    Ri_x = [ri[1]-I0[1] for ri in AllMoves]
    Ri_y = [ri[2]-I0[2] for ri in AllMoves]
    return BQ_WaveFunctionBuffer(AllMoves,psiRatios,weights,Ri_x,Ri_y)
end

function buffer_BQ_WFWeights(Walkers::Vector{<:SpiderWebWalker},ψG::AbstractGuidingFunction,Guiding_function_buffer)
    chunks = ChunkSplitters.chunks(eachindex(Walkers), n = length(Guiding_function_buffer))
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
buffer_WFWeights(O::BqOperator,Walkers,ψG,Guiding_function_buffer) = buffer_BQ_WFWeights(Walkers,ψG,Guiding_function_buffer)

measure_operator(InitialState,method::AbstractGFMCMethod,outfile,SaveConfigs,mProj,O::BqOperator,ψG::AbstractGuidingFunction,Allqs,nThreads = 2*Threads.nthreads()) = measure_operator_buffer(InitialState,method,outfile,SaveConfigs,mProj,O,ψG,Allqs,nThreads)