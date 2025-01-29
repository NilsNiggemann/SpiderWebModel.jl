"""Operator Σᵢ cos(q rᵢ) <Pᵢ+ Pᵢ†>.
"""
struct BqOperator <: AbstractOperator 
    phase::Float64
end
BqOperator() = BqOperator(0.)
BqSinOperator() = BqOperator(π/2)
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
        weights[i] = psiRatios[i]*cos((qr_i-O.phase)*0.5)^2
        wsum += weights[i]
    end
    if iszero(wsum)
        return wsum
    end
    moveidx = StatsBase.sample(StatsBase.Weights(weights,wsum))

    move1 = AllMoves[moveidx]

    applyPlaquette!(Walker.Config, move1...)
    return wsum/sqrt(length(Walker.Plaquette_positions))
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

"""Given the measurement ∑ᵢ cos²(qrᵢ/2) ⟨Bᵢ⟩, use that 2cos²(x/2)-1 = cos(x) to recover ∑ᵢ cos(qrᵢ) ⟨Bᵢ⟩."""
function recover_Bq_real(Bq_cos²::AbstractArray,Bq_cos²0::Real)
    return 2*Bq_cos² .- Bq_cos²0
end
"""Given the measurement ∑ᵢ cos²(qrᵢ/2) ⟨Bᵢ⟩, use that 2cos²(x/2)-1 = cos(x) to recover ∑ᵢ cos(qrᵢ) ⟨Bᵢ⟩. Assumes [begin] = B(q=0)."""
recover_Bq_real(Bq_cos²::AbstractArray) = recover_Bq_real(Bq_cos²,Bq_cos²[begin])

"""Given the measurement ∑ᵢ cos²((qrᵢ-π/2)/2) ⟨Bᵢ⟩, use that 2cos²(x/2)-1 = cos(x) to recover ∑ᵢ sin(qrᵢ) ⟨Bᵢ⟩."""
function recover_Bq_imag(Bq_cos²ᵩ::AbstractArray,Bq_cos²0::Real)
    return 2*Bq_cos²ᵩ .- Bq_cos²0
end

"""Given the measurement ∑ᵢ cos²((qrᵢ-π/2)/2) ⟨Bᵢ⟩, use that 2cos²(x/2)-1 = cos(x) to recover ∑ᵢ sin(qrᵢ) ⟨Bᵢ⟩. Assumes that Bq_cos²ᵩ[begin] = B_cos²ᵩ(q=0). Where B_cos²ᵩ(q=0) = 1/√2 Bq_cos²₀(0)"""
recover_Bq_imag(Bq_cos²ᵩ::AbstractArray) = recover_Bq_imag(Bq_cos²ᵩ,2Bq_cos²ᵩ[begin])

"""Recover ∑ᵢ cos(qrᵢ) ⟨Bᵢ⟩. from BCos = ∑ᵢ cos²(qrᵢ/2) ⟨Bᵢ⟩ and BSin = ∑ᵢ sin²(qrᵢ/2) ⟨Bᵢ⟩."""
function recoverBq(BCos,BSin,BCos0::Real)
    Bq_c = recover_Bq_real(BCos,BCos0)
    Bq_s = recover_Bq_real(BSin,BCos0)
    return Bq_c + 1im*Bq_s
end

"""Recover ∑ᵢ cos(qrᵢ) ⟨Bᵢ⟩. from BCos = ∑ᵢ cos²(qrᵢ/2) ⟨Bᵢ⟩ and BSin = ∑ᵢ sin²(qrᵢ/2) ⟨Bᵢ⟩. Assumes BCos[begin] = B(q=0)."""
recoverBq(BCos,BSin) = recoverBq(BCos,BSin,BCos[begin])
