struct LocalPlaquetteGuidingFunction{A<:AbstractVector} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::LocalPlaquetteGuidingFunction) = "LocalPlaquetteGuidingFunction"

function localPlaquetteGuidingFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(Float32,N)
    ψ = LocalPlaquetteGuidingFunction(params)

    get_alpha_i(ψ) .= α
    return ψ
end
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::LocalPlaquetteGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

function get_alpha_i(ψG::LocalPlaquetteGuidingFunction)
    return ψG.params
end

(ψG::LocalPlaquetteGuidingFunction)(N□::AbstractArray) = exp(guidingfunc_exponent(ψG,N□))

function guidingfunc_exponent(ψG::LocalPlaquetteGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    exponent = α' * N□
    return exponent
end

function guidingfuncRatio_exponent(ψG::LocalPlaquetteGuidingFunction,Walker::SpiderWebWalker,move,affectedPlaquettes)
    α = get_alpha_i(ψG)

    exponent = zero(eltype(α))
    n = Walker.n_x
    n´ = Walker.n_x´
    for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        exponent += α[i]*Δn
    end

    return exponent
end
