struct PlaquetteNumberGuidingFunction <: AbstractGuidingFunction
    α::Float64
end
function compute_GWF_buffer!(Buffer,ψG::PlaquetteNumberGuidingFunction,Walker::SpiderWebWalker) 
    getNPlaq!(Walker)
    return Buffer
end
get_params(ψG::PlaquetteNumberGuidingFunction) = ψG.α

(ψG::PlaquetteNumberGuidingFunction)(ΔNPlaq::Real) = exp(ψG.α*ΔNPlaq)
(ψG::PlaquetteNumberGuidingFunction)(W::SpiderWebWalker) = exp(ψG.α*sum(W.n_x))

guidingfunc_name(F::PlaquetteNumberGuidingFunction) = "PlaquetteNumberGuidingFunction"

function guidingfuncRatio_log(ψG::PlaquetteNumberGuidingFunction,Walker::SpiderWebWalker,move,AffectedPlaquetteList::T) where T

    i,j,opNum = move
    affectedPlaquettes = AffectedPlaquetteList[i,j]

    Config = get_config(Walker)
    applyPlaquette!(Config,i,j,opNum)
    getNPlaqfilled!(Walker,affectedPlaquettes)
    applyPlaquette!(Config,i,j,-opNum)
    
    α = ψG.α
    exponent = zero(α)
    n = Walker.n_x
    n´ = Walker.n_x´
    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        exponent += Δn
    end

    return exponent * α
end
# guidingfuncRatio_exponent(ψG::PlaquetteNumberGuidingFunction,Walker::SpiderWebWalker,move::Tuple,affectedPlaquettes) = _guidingfuncRatio_exponent(ψG,Walker,affectedPlaquettes)
# @inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::PlaquetteNumberGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)
