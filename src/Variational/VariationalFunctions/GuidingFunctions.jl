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

variational_parameters(P::AbstractGuidingFunction) = Dict([:params=>P.params])
variational_parameters(P::Function) = Dict([x => getproperty(P,x) for x in propertynames(P)])

@inline guidingfuncRatio(ψG::AbstractGuidingFunction,Walker::SpiderWebWalker,move,affectedPlaquettes) = exp(guidingfuncRatio_exponent(ψG,Walker,move,affectedPlaquettes))
function _guidingfuncRatio_exponent(α::AbstractVector,β::AbstractMatrix,Walker::SpiderWebWalker,affectedPlaquettes)

    exponent = zero(eltype(α))
    n = Walker.n_x
    n´ = Walker.n_x´
    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        iszero(Δn) && continue
        
        exp_i = α[i]
        LoopVectorization.@turbo for j in eachindex(n,n´)
            exp_i += β[j,i]*(n´[j] + n[j])
        end
        exponent += exp_i*Δn
    end
    return exponent
end

guidingfuncRatio_exponent(ψG::AbstractGuidingFunction,Walker::SpiderWebWalker,move) = guidingfuncRatio_exponent(ψG,eachindex(Walker.n,Walker.n´),move)

function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::AbstractGuidingFunction,Λ=0)
    (;Config,weights) = Walker
    empty!(weights)
    moves = getMoves!(Walker)
    ψx = ψG(Walker)
    for operator in moves
        i,j,opNum = operator
        applyPlaquette!(Config, i, j, opNum)
        
        ψx´ = ψG(Walker)
        
        weight = ψx´/ψx
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end

function updateWeightList_plaqs!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::AbstractGuidingFunction,Λ=0)
    (;Config,weights) = Walker
    empty!(weights)
    moves = getMoves!(Walker)
    n_x = getNPlaq!(Walker)
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        n_x´ = getNPlaqfilled!(Walker,indices)
        
        weight = guidingfuncRatio(ψG,Walker,operator,indices)
        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end
    if Λ != 0
        push!(moves, (0,0,0))
        push!(weights,Λ)
    end
    return weights
end


function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::T,Λ=0) where T
    (;Config,weights,moves,n_x,n_x´) = Walker
    isempty(weights) && isempty(moves) || return weights #return if weights are already computed
    
    getMoves!(Walker)
    # getNPlaq!(Walker)
    
    psix_inv = 1/ψG(get_config(Walker))
    
    for operator in moves
        i,j,opNum = operator
        indices = AffectedPlaquetteList[i,j]
        applyPlaquette!(Config, i, j, opNum)
        
        # getNPlaq!(Walker,indices)
        
        # N□ = getNPlaq_difference(n_x,n_x´,indices) 
        weight = ψG(get_config(Walker)) * psix_inv

        push!(weights,weight)
        applyPlaquette!(Config, i, j, -opNum)
    end

    
    if Λ != 0
        push!(moves, DIAGONAL_MOVE_ID)
        push!(weights,Λ)
    end
    return weights
end