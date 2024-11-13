abstract type AbstractGuidingFunction end

"""Contains guiding wave function and information about symmetry equivalent parameters"""
struct SymmetryReducedWaveFunction{F<:AbstractGuidingFunction} <: AbstractGuidingFunction
    psi::F
    indicesMapping::Vector{Int}
    uniqueInds::Vector{Int}
end

(ψG::AbstractGuidingFunction)(x::AbstractWalker) = ψG(get_config(x))

guidingfunc_name(F::Function) = string(typeof(F))
guidingfunc_name(F::AbstractGuidingFunction) = string(typeof(F))

get_params(P::AbstractGuidingFunction) = P.params
get_params(P::SymmetryReducedWaveFunction) = get_params(P.psi)
variational_parameters(P::AbstractGuidingFunction) = Dict([:params=>P.params])
variational_parameters(P::Function) = Dict([x => getproperty(P,x) for x in propertynames(P)])

"""Allocates a buffer for the guiding wave function to be used for efficient computation. Can be anything depending on the model"""
allocate_GWF_buffer(ψG::AbstractGuidingFunction,S::Any) = precomputeAffectedPlaquettes(S)

function allocate_GWF_buffers_threads(ψG::AbstractGuidingFunction,InitialState)
    fetch.([Threads.@spawn allocate_GWF_buffer(ψG, InitialState) for _ in 1:Threads.nthreads()])
end

"""precomputes the guiding wave function for the given configuration. Can be overloaded for more complex wavefunctions in which case a Buffer is filled. In this case the buffer needs to be returned"""
fill_GWF_buffer!(Buff,ψG::AbstractGuidingFunction,Walker) = ψG(Walker)

# allocate_GWF_buffer(ψG::AbstractGuidingFunction,S::StencilSpinConfig) = nothing

function updateWeightList!(Walker::SpiderWebWalker,Buffer,ψG::AbstractGuidingFunction)
    weights = Walker.weights
    empty!(weights)
    moves = getMoves!(Walker)
    newBuff = fill_GWF_buffer!(Buffer,ψG,Walker)

    for operator in moves
        weight = guidingfuncRatio(ψG,Walker,operator,newBuff)
        push!(weights,weight)
    end

    return weights
end

function updateWeightList_naive!(Walker::SpiderWebWalker,Buffer,ψG::AbstractGuidingFunction)
    weights = Walker.weights
    empty!(weights)
    moves = getMoves!(Walker)
    for operator in moves
        psix = ψG(Walker)
        i,j,opsign = operator
        applyPlaquette!(get_config(Walker),i,j,opsign)
        psix´ = ψG(Walker)
        applyPlaquette!(get_config(Walker),i,j,-opsign)
        weight = psix´/psix

        push!(weights,weight)
    end

    return weights
end

function updateWeightList!(Walker::SpiderWebWalker,Buffer,ψG::AbstractGuidingFunction,Λ::Real)
    res = updateWeightList!(Walker,Buffer,ψG)
    if Λ != 0
        push!(Walker.moves, DIAGONAL_MOVE_ID)
        push!(Walker.weights,Λ)
    end
    return res
end

function guidingfuncRatio(ψG::AbstractGuidingFunction,Walker,operator,ψx::Number)
    i,j,opNum = operator
    Config = get_config(Walker)
    applyPlaquette!(Config, i, j, opNum)
    ψx´ = ψG(Walker)
    applyPlaquette!(Config, i, j, -opNum)
    return ψx´/ψx
end

function _getParamsTypeAndIndex(A::RecursiveArrayTools.ArrayPartition,i)
    @boundscheck checkbounds(A, i)
    @inbounds for j in 1:length(A.x)
        i -= length(A.x[j])
        if i <= 0
            return (j,length(A.x[j]) + i)
        end
    end
    throw(BoundsError(A, i))
end