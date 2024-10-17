struct RKFunction <: AbstractGuidingFunction end

guidingfunc_name(F::RKFunction) = "RKFunction"
variational_parameters(::RKFunction) = Dict{Symbol,Float64}()
guidingfuncRatio_exponent(ψG::RKFunction,Walker::SpiderWebWalker,affectedPlaquettes) = 0.
function updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::RKFunction,Λ=0)
    (;Config,weights,moves,n_x,n_x´) = Walker
    isempty(weights) && isempty(moves) || return weights #return if weights are already computed
    
    getMoves!(Walker)
    
    resize!(weights,length(moves))
    weights .= 1.

    if Λ != 0
        push!(moves, DIAGONAL_MOVE_ID)
        push!(weights,Λ)
    end
    return weights
end