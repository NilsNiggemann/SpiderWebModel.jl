struct FullVariationalGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
end
guidingfunc_name(F::FullVariationalGuidingFunction) = "FullVariationalGuidingFunction"

function fullVariationalFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    params = zeros(N,N+1)
    ψ = FullVariationalGuidingFunction(params)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin]
end

function get_beta_ij(ψG::FullVariationalGuidingFunction)
    return @view ψG.params[:,begin+1:end]
end

function (ψG::FullVariationalGuidingFunction)(N□::AbstractArray) 
    return exp(guidingfunc_exponent(ψG,N□))
end

function guidingfunc_exponent(ψG::FullVariationalGuidingFunction,N□::AbstractArray) 
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    exponent = α' * N□ + dot(N□,β,N□)
    return exponent
end

@inline guidingfuncRatio(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,move) = exp(guidingfuncRatio_exponent(ψG,n,n´,move))

function guidingfuncRatio_exponent(ψG::FullVariationalGuidingFunction,Walker::SpiderWebWalker,move,affectedPlaquettes)
    α = get_alpha_i(ψG)
    β = get_beta_ij(ψG)
    _guidingfuncRatio_exponent(α,β,Walker,affectedPlaquettes)
end

@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::FullVariationalGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)
