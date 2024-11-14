
struct JastrowFunction{T<:Number} <: AbstractGuidingFunction
    α::Vector{T}
    m_i::Vector{T}
    v_ij::Matrix{T}
end
function JastrowFunction(N::Int,NPlaq,Type = Float32)
    α = zeros(Type,NPlaq)
    m_i = zeros(Type,N)
    v_ij = zeros(Type,N,N)
    return JastrowFunction(α,m_i,v_ij)
end
JastrowFunction(S::StencilSpinConfig,Type = Float32) = JastrowFunction(length(S),length(collect(plaquetteIterator(S))),Type)

get_alpha_i(ψG::JastrowFunction) = ψG.α
get_m_i(ψG::JastrowFunction) = ψG.m_i
get_v_ij(ψG::JastrowFunction) = ψG.v_ij

function get_params(ψG::JastrowFunction)
    return RecursiveArrayTools.ArrayPartition(ψG.α,ψG.m_i,ψG.v_ij)
end

function (ψG::JastrowFunction)(x::StencilSpinConfig)
    n = getNPlaq(x)
    _evaluate_jastrow(ψG,x,n)
end

guidingfunc_name(F::JastrowFunction) = "JastrowFunction"
(ψG::JastrowFunction)(W::SpiderWebWalker) = _evaluate_jastrow(ψG,get_config(W),W.n_x)

function _evaluate_jastrow(ψG,x::AbstractMatrix,n)
    α = get_alpha_i(ψG)
    m = get_m_i(ψG)
    W = get_v_ij(ψG)

    exp_α = zero(eltype(α))
    for i in eachindex(n)
        exp_α += α[i] * n[i]
    end
    exp_m = zero(eltype(α))
    for i in eachindex(IndexLinear(),x)
        exp_m += m[i] * x[i]
    end
    x_lin = reshape(x, length(x))
    exp_W = dot(x_lin,W,x_lin)
    return exp(exp_α + exp_m + exp_W)
end

struct Jastrow_GWF_Buffer_3{T1<:AbstractVector,D,L}
    h_i::T1
    prefac_moves::D
    AffectedPlaquetteList::L
end

function allocate_GWF_buffer(ψG::JastrowFunction,S::AbstractMatrix)
    h_i = zeros(length(S))
    prefac_moves = _precompute_prefac_moves(ψG,S)
    AffectedPlaquetteList = precomputeAffectedPlaquettes(S)
    return Jastrow_GWF_Buffer_3(h_i,prefac_moves,AffectedPlaquetteList)
end


function _precompute_prefac_moves(ψG::JastrowFunction,Conf::StencilSpinConfig)
    AllPlaqs = collect(plaquetteIterator(Conf))
    vij = get_v_ij(ψG)
    AllWeights = Dict((i,j) => _precompute_jastrow_weight(vij,(i,j),Conf) for (i,j) in AllPlaqs)
end

function _precompute_jastrow_weight(vij,move,Conf::StencilSpinConfig)
    LI = LinearIndices(Conf)
    i,j = move

    sites = safe_parent_indices(parent(Conf), (i,j))

    exponent = zero(eltype(vij))

    for k in eachindex(P1_STENCIL)
        F_k = P1_STENCIL[k]#*opSign
        a_k = LI[CartesianIndex(sites[k])]

        for k´ in eachindex(P1_STENCIL)
            F_k´  = P1_STENCIL[k´]#*opSign (will cancel out)
            a_k´ =  LI[CartesianIndex(sites[k´])]
            exponent += 0.5*vij[a_k,a_k´]*F_k*F_k´
        end
    end
    return exp(exponent)
end


function fill_GWF_buffer!(Buffer,ψG::JastrowFunction,Walker::SpiderWebWalker) 
    getNPlaq!(Walker)
    x = get_config(Walker)
    v = get_v_ij(ψG)
    for j in axes(v,2)
        for i in axes(v,1)
            Buffer.h_i[j] = x[i]*v[i,j]
        end
    end
    return Buffer
end


function guidingfuncRatio(ψG::JastrowFunction,Walker::SpiderWebWalker,move::Tuple,Buffer::Jastrow_GWF_Buffer_3)
    α = get_alpha_i(ψG)
    m = get_m_i(ψG)

    (;h_i,prefac_moves,AffectedPlaquetteList) = Buffer
    
    x = get_config(Walker)

    i,j,opSign = move
    affectedPlaquettes = AffectedPlaquetteList[i,j]
    prefac = prefac_moves[(i,j)]

    sites = safe_parent_indices(parent(x), (i, j))

    LI = LinearIndices(x)

    applyPlaquette!(x, i, j, opSign)
    getNPlaqfilled!(Walker,affectedPlaquettes)
    applyPlaquette!(x, i, j, -opSign)

    exp_α = zero(eltype(get_params(ψG)))
    exp_h = zero(eltype(get_params(ψG)))
    exp_m = zero(eltype(get_params(ψG)))

    n = Walker.n_x
    n´ = Walker.n_x´

    @inbounds @simd for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        exp_α += α[i] * Δn
    end

    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign
        I = LI[i,j]
        exp_h += h_i[I]*s
        exp_m += m[I]*x[I]
    end

    return exp(exp_α + exp_h + exp_m)
end
function getOx_k(ψG::JastrowFunction,W::SpiderWebWalker,k)
    par = get_params(ψG)
    x = get_config(W)
    n = W.n_x
    paramtype,k = _getParamsTypeAndIndex(par,k)

    if paramtype == 1 
        return n[k]
    end

    if paramtype == 2
        return x[k]
    end

    vij = get_v_ij(ψG)
    
    if paramtype == 3
        i,j = Tuple(CartesianIndices(vij)[k])
        return x[i] * x[j]
    end
    error("Invalid parameter type")
end
