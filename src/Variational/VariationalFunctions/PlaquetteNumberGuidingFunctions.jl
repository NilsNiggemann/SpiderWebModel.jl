struct PlaquetteNumberGuidingFunction <: AbstractGuidingFunction
    α::Float64
end
(ψG::PlaquetteNumberGuidingFunction)(ΔNPlaq::Real) = exp(ψG.α*ΔNPlaq)
guidingfunc_name(F::PlaquetteNumberGuidingFunction) = "PlaquetteNumberGuidingFunction"
variational_parameters(P::PlaquetteNumberGuidingFunction) = Dict([:alpha=>P.α])
function guidingfuncRatio_exponent(ψG::PlaquetteNumberGuidingFunction,Walker::SpiderWebWalker,affectedPlaquettes)
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
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::PlaquetteNumberGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)
