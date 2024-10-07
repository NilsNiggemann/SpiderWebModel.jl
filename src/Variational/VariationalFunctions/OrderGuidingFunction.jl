struct OrderGuidingFunction{A<:AbstractArray} <: AbstractGuidingFunction
    params::A
    N::Int
    Nsites::Int
end
guidingfunc_name(F::OrderGuidingFunction) = "OrderGuidingFunction"
function orderGuidingFunction(State,α::Real=0.1)
    plaqs = collect(plaquetteIterator(State))
    N = length(plaqs)
    Nsites = length(State)
    params = zeros(N + Nsites)
    ψ = OrderGuidingFunction(params,N,Nsites)
    get_alpha_i(ψ) .= α
    return ψ
end

function get_alpha_i(ψG::OrderGuidingFunction)
    @view ψG.params[begin:begin+ψG.N-1]
end

function get_m_i(ψG::OrderGuidingFunction)
    startIndex = firstindex(ψG.params) + ψG.N
    m = @view ψG.params[startIndex:startIndex + ψG.Nsites-1]
    return m
end

function guidingfuncRatio_exponent(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,move::Tuple,affectedPlaquettes)
    α = get_alpha_i(ψG)
    m = get_m_i(ψG)

    Mat = parent(get_config(Walker))

    exp_α = zero(eltype(ψG.params))
    n = Walker.n_x
    n´ = Walker.n_x´
    for i in affectedPlaquettes
        Δn = n´[i] - n[i]
        exp_α += α[i] * Δn
    end

    i,j,opSign = move
    
    sites = safe_parent_indices(Mat, (i, j))
    LI = LinearIndices(Mat)

    exp_m = zero(exp_α)
    for (ij,s) in zip(sites,P1_STENCIL)
        i,j = ij
        I = LI[i,j]
        exp_m += m[I]*s
    end
    return exp_α + opSign*exp_m
end
@inline updateWeightList!(Walker::SpiderWebWalker,AffectedPlaquetteList,ψG::OrderGuidingFunction,Λ=0) = updateWeightList_plaqs!(Walker,AffectedPlaquetteList,ψG,Λ)

function getOx_k(ψG::OrderGuidingFunction,Walker::SpiderWebWalker,k)
    α = get_alpha_i(ψG)

    if k in eachindex(α)
        return Walker.n_x[k]
    end
    # m = get_m_i(ψG)

    Ok = get_config(Walker)[k-lastindex(α)]

    return Ok
end
