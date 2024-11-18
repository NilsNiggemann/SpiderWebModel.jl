
struct JastrowFunction{T<:Number} <: AbstractGuidingFunction
    α::Vector{T}
    m_i::Vector{T}
    v_ij::Matrix{T}
end
function JastrowFunction(N::Int,NPlaq,Type = Float64)
    α = zeros(Type,NPlaq)
    m_i = zeros(Type,N)
    v_ij = zeros(Type,N,N)
    return JastrowFunction(α,m_i,v_ij)
end
JastrowFunction(S::StencilSpinConfig,Type = Float64) = JastrowFunction(length(S),length(collect(plaquetteIterator(S))),Type)

get_alpha_i(ψG::JastrowFunction) = ψG.α
get_m_i(ψG::JastrowFunction) = ψG.m_i
get_v_ij(ψG::JastrowFunction) = ψG.v_ij

function get_params(ψG::JastrowFunction)
    return RecursiveArrayTools.ArrayPartition(ψG.α,ψG.m_i,get_v_ij(ψG))
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
    v = get_v_ij(ψG)

    exp_α = zero(eltype(α))
    for i in eachindex(n)
        exp_α += α[i] * n[i]
    end
    exp_m = zero(eltype(α))
    for i in eachindex(IndexLinear(),x)
        exp_m += m[i] * x[i]
    end
    x_lin = reshape(x, length(x))
    exp_v = 0.5*dot(x_lin,v,x_lin)

    return exp(exp_α + exp_m + exp_v)
end

struct Jastrow_GWF_Buffer{T<:Number,D,L}
    x_i::Vector{T}
    h_i::Vector{T}
    prefac_moves::D
    AffectedPlaquetteList::L
end

function allocate_GWF_buffer(ψG::JastrowFunction{T},S::AbstractMatrix) where T
    x_i = zeros(T,length(S))
    h_i = zeros(T,length(S))
    
    prefac_moves = _precompute_prefac_moves(ψG,S)
    AffectedPlaquetteList = precomputeAffectedPlaquettes(S)
    return Jastrow_GWF_Buffer(x_i,h_i,prefac_moves,AffectedPlaquetteList)
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
    return exponent
end


function premove_update_GWF_buffer!(Buffer::Jastrow_GWF_Buffer,ψG::JastrowFunction,Walker::SpiderWebWalker) 
    getNPlaq!(Walker)
    return Buffer
end

function compute_GWF_buffer!(Buffer::Jastrow_GWF_Buffer,ψG::JastrowFunction,Walker::SpiderWebWalker)

    getNPlaq!(Walker)
    x = Buffer.x_i

    x .= reshape(get_config(Walker),length(x))

    v = get_v_ij(ψG)
    h = Buffer.h_i
    # fill!(h,0.0)

    for i in axes(v,1)
        hi = 0.
        LoopVectorization.@turbo for j in axes(v,2)
            hi += x[j]*v[j,i]
        end
        h[i] = hi
    end
    # mul!(h,v,x) 
    return Buffer
end

function post_move_update_GWF_buffer!(Buffer::Jastrow_GWF_Buffer,ψG::JastrowFunction,Walker::SpiderWebWalker,move::Tuple)
    Config = get_config(Walker)

    i,j,opSign = move

    affected_sites = safe_parent_indices(parent(Config), (i, j))

    v = get_v_ij(ψG)
    h = Buffer.h_i

    LI = LinearIndices(Config)
    
    for (index,site) in enumerate(affected_sites)

        i = LI[CartesianIndex(site)]
        s = P1_STENCIL[index]*opSign

        for j in eachindex(h)
            h[j] += v[j,i]*s
        end
    end
    return Buffer
end

function guidingfuncRatio(ψG::JastrowFunction,Walker::SpiderWebWalker,move::Tuple,Buffer::Jastrow_GWF_Buffer)
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

    LoopVectorization.@turbo for ind in eachindex(affectedPlaquettes)
        i = affectedPlaquettes[ind]
        Δn = n´[i] - n[i]
        exp_α += α[i] * Δn
    end

    @inbounds @simd for idx in eachindex(sites)
        i,j = sites[idx]
        s = P1_STENCIL[idx]*opSign
        I = LI[i,j]
        exp_h += h_i[I]*s
        exp_m += m[I]*s
    end

    return exp(exp_α + exp_h + exp_m + prefac)
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
        return 0.5*x[i] * x[j]
    end
    error("Invalid parameter type")
end

function generate_equivalent(type,k,T::TranslationalSymmetry,ψ::JastrowFunction,S::AbstractMatrix) 
    
    if type == 1
        AllPlaqs = collect(plaquetteIterator(S))

        site = AllPlaqs[k]
        sites = generate_equivalent_sites(site,T,S)
        
        return plaquette_to_index.(sites,Ref(AllPlaqs))
    elseif type == 2
        site = index_to_site(k,S)
        sites = generate_equivalent_sites(site,T,S)
        return site_to_index.(sites,Ref(S))
    elseif type == 3
        siteI,siteJ = index_to_site_pair(k,S,ψ.v_ij)
        sites = generate_equivalent_site_pairs(siteI,siteJ,T,S)
        i_sites = getindex.(sites,1)
        j_sites = getindex.(sites,2)
        return site_pair_to_index.(i_sites,j_sites,Ref(S),Ref(ψ.v_ij))
    else
        error("Invalid type")
    end

end

struct ExchangeSymmetry <: AbstractSymop end

function generate_equivalent(type,k,T::ExchangeSymmetry,ψ::JastrowFunction,S::AbstractMatrix) 
    if type == 3
        i,j = index_to_site_pair(k,S,ψ.v_ij)
        return [k,site_pair_to_index(j,i,S,ψ.v_ij)]
    else
        return [k]
    end
end